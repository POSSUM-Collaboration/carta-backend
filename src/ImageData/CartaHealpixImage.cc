/* This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-backend
   Copyright 2018- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
   Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
   SPDX-License-Identifier: GPL-3.0-or-later
*/

//# CartaHealpixImage.cc : HEALPix map presented as a casacore::ImageInterface via on-the-fly reprojection

#include "CartaHealpixImage.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <list>
#include <mutex>
#include <unordered_map>

#include <fcntl.h>
#include <fitsio.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

#include <casacore/casa/Arrays/Slice.h>
#include <casacore/casa/Exceptions/Error.h>
#include <casacore/casa/Quanta/MVDirection.h>
#include <casacore/casa/Quanta/Quantum.h>
#include <casacore/casa/Quanta/UnitMap.h>
#include <casacore/casa/Quanta/UnitVal.h>
#include <casacore/coordinates/Coordinates/LinearCoordinate.h>
#include <casacore/coordinates/Coordinates/Projection.h>
#include <casacore/coordinates/Coordinates/SpectralCoordinate.h>
#include <casacore/coordinates/Coordinates/StokesCoordinate.h>
#include <casacore/images/Images/ImageInfo.h>
#include <casacore/measures/Measures/MCDirection.h>
#include <casacore/measures/Measures/MFrequency.h>
#include <casacore/measures/Measures/MeasConvert.h>
#include <casacore/measures/Measures/Stokes.h>
#include <casacore/scimath/Mathematics/GaussianBeam.h>

#include <healpix_base.h>
#include <pointing.h>

#include "Logger/Logger.h"
#include "Util/Nan.h"

using namespace carta;

namespace {

const double PI = M_PI;
const double DEG2RAD = M_PI / 180.0;
const unsigned int CURSOR_SIZE = 512;
// Most pixels per side of a grid: the tile coordinates of the image viewer index 4096 tiles of 256 pixels per axis
const int MAX_GRID_SIDE = 1 << 20;

std::string UpperTrim(std::string s) {
    s.erase(0, s.find_first_not_of(" \t'"));
    s.erase(s.find_last_not_of(" \t'") + 1);
    std::transform(s.begin(), s.end(), s.begin(), ::toupper);
    return s;
}

bool ReadStringKey(fitsfile* fptr, const char* key, std::string& value) {
    int status(0);
    char buffer[FLEN_VALUE];
    fits_read_key(fptr, TSTRING, key, buffer, nullptr, &status);
    if (status) {
        return false;
    }
    value = UpperTrim(buffer);
    return true;
}

bool ReadLongKey(fitsfile* fptr, const char* key, long& value) {
    int status(0);
    fits_read_key(fptr, TLONG, key, &value, nullptr, &status);
    return status == 0;
}

bool ReadDoubleKey(fitsfile* fptr, const char* key, double& value) {
    int status(0);
    fits_read_key(fptr, TDOUBLE, key, &value, nullptr, &status);
    return status == 0;
}

// Map a HEALPix table column name to a Stokes type, or Stokes::Undefined
casacore::Stokes::StokesTypes ColumnStokes(const std::string& column) {
    std::string name = UpperTrim(column);
    if (name == "I" || name == "TEMPERATURE" || name == "I_STOKES" || name == "INTENSITY" || name == "T") {
        return casacore::Stokes::I;
    } else if (name == "Q" || name == "Q_POLARISATION" || name == "Q_POLARIZATION" || name == "Q_STOKES") {
        return casacore::Stokes::Q;
    } else if (name == "U" || name == "U_POLARISATION" || name == "U_POLARIZATION" || name == "U_STOKES") {
        return casacore::Stokes::U;
    } else if (name == "V" || name == "V_POLARISATION" || name == "V_POLARIZATION" || name == "V_STOKES") {
        return casacore::Stokes::V;
    }
    return casacore::Stokes::Undefined;
}

casacore::MDirection::Types FrameFromCode(const std::string& code, bool& ok) {
    ok = true;
    std::string c = UpperTrim(code);
    if (c.empty()) {
        ok = false;
        return casacore::MDirection::ICRS;
    }
    if (c[0] == 'G') {
        return casacore::MDirection::GALACTIC;
    } else if (c[0] == 'E') {
        return casacore::MDirection::ECLIPTIC;
    } else if (c[0] == 'C' || c[0] == 'Q') {
        return casacore::MDirection::ICRS;
    } else if (c[0] == 'S') {
        return casacore::MDirection::SUPERGAL;
    }
    // Try a casacore name (J2000, FK5, B1950, ...)
    casacore::MDirection::Types type;
    if (casacore::MDirection::getType(type, c)) {
        return type;
    }
    ok = false;
    return casacore::MDirection::ICRS;
}

bool IsBadValue(float v, float bad) {
    return !std::isfinite(v) || v == bad || v < -1.0e30f;
}

} // namespace

// ---------------------------------------------------------------------------
// Map sources: where the values of each map come from

namespace carta {

/**
 * @brief Values of the maps of a HEALPix HDU, indexed by map (table column or image row) and HEALPix pixel.
 * ReadRange reads a contiguous run of pixels of one map; ReadPixels reads an ascending list of pixels, as few
 * reads as the source allows (a spectral profile reads the same pixels from every map, so this is the hot path
 * when the file is on a slow or network file system). Both must be safe to call from several threads at once.
 * Values are NaN where unseen or out of range.
 */
class HealpixMapSource {
public:
    virtual ~HealpixMapSource() = default;
    virtual size_t NumMaps() const = 0;
    virtual size_t Npix() const = 0;
    virtual void ReadRange(size_t map, int64_t first, size_t count, float* out) = 0;
    virtual void ReadPixels(size_t map, const int64_t* pixels, size_t count, float* out) = 0;
};

namespace {

// Maps held in memory (binary-table columns)
class InMemoryMapSource : public HealpixMapSource {
public:
    explicit InMemoryMapSource(std::vector<std::shared_ptr<std::vector<float>>> maps) : _maps(std::move(maps)) {}
    size_t NumMaps() const override {
        return _maps.size();
    }
    size_t Npix() const override {
        return _maps.empty() ? 0 : _maps[0]->size();
    }
    void ReadRange(size_t map, int64_t first, size_t count, float* out) override {
        if (map >= _maps.size()) {
            std::fill(out, out + count, FLOAT_NAN);
            return;
        }
        const auto& values = *_maps[map];
        for (size_t i = 0; i < count; ++i) {
            int64_t pix = first + i;
            out[i] = (pix >= 0 && pix < static_cast<int64_t>(values.size())) ? values[pix] : FLOAT_NAN;
        }
    }
    void ReadPixels(size_t map, const int64_t* pixels, size_t count, float* out) override {
        if (map >= _maps.size()) {
            std::fill(out, out + count, FLOAT_NAN);
            return;
        }
        const auto& values = *_maps[map];
        for (size_t i = 0; i < count; ++i) {
            int64_t pix = pixels[i];
            out[i] = (pix >= 0 && pix < static_cast<int64_t>(values.size())) ? values[pix] : FLOAT_NAN;
        }
    }

private:
    std::vector<std::shared_ptr<std::vector<float>>> _maps;
};

// Image-format maps read directly from an uncompressed FITS file (BITPIX -32 or -64, big-endian) with pread:
// explicit reads of the byte ranges needed, rather than page faults on a memory mapping, so that a file on a
// network or virtualised file system (a Docker or container volume) is read in a few requests per map
class MappedImageMapSource : public HealpixMapSource {
public:
    MappedImageMapSource(
        const std::string& filename, long long data_offset, int bitpix, long npix, long num_maps, double bscale, double bzero, float bad)
        : _npix(npix),
          _num_maps(num_maps),
          _bytes(std::abs(bitpix) / 8),
          _data_offset(data_offset),
          _bscale(bscale),
          _bzero(bzero),
          _bad(bad) {
        _fd = open(filename.c_str(), O_RDONLY);
        if (_fd < 0) {
            throw(casacore::AipsError("Cannot open HEALPix file for reading."));
        }
        struct stat st;
        if (fstat(_fd, &st) != 0 || st.st_size < data_offset + static_cast<long long>(_bytes) * npix * num_maps) {
            close(_fd);
            throw(casacore::AipsError("HEALPix file is shorter than its header describes."));
        }
        const uint16_t one(1);
        _swap = (*reinterpret_cast<const uint8_t*>(&one) == 1); // little-endian host, FITS data is big-endian
    }
    ~MappedImageMapSource() override {
        if (_fd >= 0) {
            close(_fd);
        }
    }
    size_t NumMaps() const override {
        return _num_maps;
    }
    size_t Npix() const override {
        return _npix;
    }
    void ReadRange(size_t map, int64_t first, size_t count, float* out) override {
        if (map >= static_cast<size_t>(_num_maps) || first < 0 || first + static_cast<int64_t>(count) > _npix) {
            std::fill(out, out + count, FLOAT_NAN);
            return;
        }
        std::vector<char> raw(count * _bytes);
        if (!ReadBytes(map, first, raw)) {
            std::fill(out, out + count, FLOAT_NAN);
            return;
        }
        for (size_t i = 0; i < count; ++i) {
            out[i] = Decode(raw.data() + i * _bytes);
        }
    }
    void ReadPixels(size_t map, const int64_t* pixels, size_t count, float* out) override {
        // Pixels are ascending: read them in runs, starting a new run where the gap would cost more than a
        // small read or the run's buffer would get large; each run is one pread of a contiguous byte range
        const int64_t max_gap(2048);
        const int64_t max_run = (4 << 20) / _bytes;
        std::vector<char> raw;
        size_t start(0);
        while (start < count) {
            size_t end(start + 1);
            while (end < count && pixels[end] - pixels[end - 1] <= max_gap && pixels[end] - pixels[start] < max_run) {
                ++end;
            }
            if (map >= static_cast<size_t>(_num_maps) || pixels[start] < 0 || pixels[end - 1] >= _npix) {
                for (size_t i = start; i < end; ++i) {
                    out[i] =
                        (pixels[i] >= 0 && pixels[i] < _npix && map < static_cast<size_t>(_num_maps)) ? ReadOne(map, pixels[i]) : FLOAT_NAN;
                }
            } else {
                const size_t run = pixels[end - 1] - pixels[start] + 1;
                raw.resize(run * _bytes);
                if (ReadBytes(map, pixels[start], raw)) {
                    for (size_t i = start; i < end; ++i) {
                        out[i] = Decode(raw.data() + (pixels[i] - pixels[start]) * _bytes);
                    }
                } else {
                    std::fill(out + start, out + end, FLOAT_NAN);
                }
            }
            start = end;
        }
    }

private:
    bool ReadBytes(size_t map, int64_t first, std::vector<char>& raw) const {
        off_t offset = _data_offset + (static_cast<int64_t>(map) * _npix + first) * _bytes;
        size_t done(0);
        while (done < raw.size()) {
            ssize_t n = pread(_fd, raw.data() + done, raw.size() - done, offset + done);
            if (n <= 0) {
                return false;
            }
            done += n;
        }
        return true;
    }
    float ReadOne(size_t map, int64_t pix) const {
        std::vector<char> raw(_bytes);
        return ReadBytes(map, pix, raw) ? Decode(raw.data()) : FLOAT_NAN;
    }
    float Decode(const char* p) const {
        double v;
        if (_bytes == 4) {
            uint32_t u;
            std::memcpy(&u, p, 4);
            if (_swap) {
                u = __builtin_bswap32(u);
            }
            float f;
            std::memcpy(&f, &u, 4);
            v = f;
        } else {
            uint64_t u;
            std::memcpy(&u, p, 8);
            if (_swap) {
                u = __builtin_bswap64(u);
            }
            double d;
            std::memcpy(&d, &u, 8);
            v = d;
        }
        float value = static_cast<float>(v * _bscale + _bzero);
        return IsBadValue(value, _bad) ? FLOAT_NAN : value;
    }

    long _npix;
    long _num_maps;
    int _bytes;
    long long _data_offset;
    double _bscale, _bzero;
    float _bad;
    int _fd = -1;
    bool _swap = true;
};

// Image-format maps read a whole map at a time through cfitsio (compressed files, integer data), with a small cache
class CfitsioRowMapSource : public HealpixMapSource {
public:
    CfitsioRowMapSource(const std::string& filename, unsigned int hdu, long npix, long num_maps, float bad)
        : _filename(filename), _hdu(hdu), _npix(npix), _num_maps(num_maps), _bad(bad) {}
    ~CfitsioRowMapSource() override {
        if (_fptr) {
            int status(0);
            fits_close_file(_fptr, &status);
        }
    }
    size_t NumMaps() const override {
        return _num_maps;
    }
    size_t Npix() const override {
        return _npix;
    }
    void ReadRange(size_t map, int64_t first, size_t count, float* out) override {
        auto values = FullMap(map);
        for (size_t i = 0; i < count; ++i) {
            int64_t pix = first + i;
            out[i] = (values && pix >= 0 && pix < static_cast<int64_t>(values->size())) ? (*values)[pix] : FLOAT_NAN;
        }
    }
    void ReadPixels(size_t map, const int64_t* pixels, size_t count, float* out) override {
        auto values = FullMap(map);
        for (size_t i = 0; i < count; ++i) {
            int64_t pix = pixels[i];
            out[i] = (values && pix >= 0 && pix < static_cast<int64_t>(values->size())) ? (*values)[pix] : FLOAT_NAN;
        }
    }

private:
    std::shared_ptr<const std::vector<float>> FullMap(size_t map) {
        if (map >= static_cast<size_t>(_num_maps)) {
            return nullptr;
        }
        std::lock_guard<std::mutex> guard(_mutex);
        auto found = _cache.find(map);
        if (found != _cache.end()) {
            _order.remove(map);
            _order.push_front(map);
            return found->second;
        }
        if (!_fptr) {
            int status(0);
            fits_open_file(&_fptr, _filename.c_str(), READONLY, &status);
            int hdu_type(0);
            fits_movabs_hdu(_fptr, _hdu + 1, &hdu_type, &status);
            if (status) {
                throw(casacore::AipsError("Cannot open HEALPix image HDU."));
            }
        }
        auto values = std::make_shared<std::vector<float>>(_npix);
        int status(0), anynul(0);
        float nulval = FLOAT_NAN;
        fits_read_img(_fptr, TFLOAT, static_cast<LONGLONG>(map) * _npix + 1, _npix, &nulval, values->data(), &anynul, &status);
        if (status) {
            throw(casacore::AipsError("Error reading HEALPix map " + std::to_string(map) + "."));
        }
        for (auto& v : *values) {
            if (IsBadValue(v, _bad)) {
                v = FLOAT_NAN;
            }
        }
        _cache[map] = values;
        _order.push_front(map);
        while (_order.size() > MAX_CACHED_MAPS) {
            _cache.erase(_order.back());
            _order.pop_back();
        }
        return values;
    }

    static constexpr size_t MAX_CACHED_MAPS = 8;
    std::string _filename;
    unsigned int _hdu;
    long _npix;
    long _num_maps;
    float _bad;
    fitsfile* _fptr = nullptr;
    std::mutex _mutex;
    std::unordered_map<size_t, std::shared_ptr<std::vector<float>>> _cache;
    std::list<size_t> _order; // most recently used first
};

bool IsCompressedFilename(const std::string& filename) {
    for (const char* suffix : {".gz", ".GZ", ".bz2", ".Z", ".zip", ".fz"}) {
        size_t n = std::strlen(suffix);
        if (filename.size() > n && filename.compare(filename.size() - n, n, suffix) == 0) {
            return true;
        }
    }
    return false;
}

} // namespace
} // namespace carta

// ---------------------------------------------------------------------------
// HealpixMapInfo

long HealpixMapInfo::npix() const {
    return 12 * nside * nside;
}

double HealpixMapInfo::resolution_deg() const {
    return std::sqrt(4.0 * PI / npix()) / DEG2RAD;
}

casacore::MDirection::Types HealpixMapInfo::frame() const {
    bool ok;
    auto type = FrameFromCode(coordsys, ok);
    return ok ? type : casacore::MDirection::ICRS;
}

// ---------------------------------------------------------------------------
// Static metadata helpers

bool CartaHealpixImage::IsHealpixHdu(const std::string& filename, unsigned int hdu, HealpixMapInfo* info) {
    std::string error;
    HealpixMapInfo local_info;
    bool is_healpix = ReadMapInfo(filename, hdu, info ? *info : local_info, error);
    return is_healpix;
}

bool CartaHealpixImage::IsHealpixFile(const std::string& filename, unsigned int* first_hdu) {
    // HEALPix maps are usually binary tables in an extension, or an image HDU with one map per row; check the first few HDUs.
    fitsfile* fptr(nullptr);
    int status(0);
    fits_open_file(&fptr, filename.c_str(), READONLY, &status);
    if (status) {
        return false;
    }
    int num_hdus(0);
    fits_get_num_hdus(fptr, &num_hdus, &status);
    fits_close_file(fptr, &status);

    for (int hdu = 0; hdu < std::min(num_hdus, 4); ++hdu) {
        if (IsHealpixHdu(filename, hdu)) {
            if (first_hdu) {
                *first_hdu = hdu;
            }
            return true;
        }
    }
    return false;
}

bool CartaHealpixImage::ReadMapInfo(const std::string& filename, unsigned int hdu, HealpixMapInfo& info, std::string& error) {
    fitsfile* fptr(nullptr);
    int status(0);
    fits_open_file(&fptr, filename.c_str(), READONLY, &status);
    if (status) {
        error = "Cannot open FITS file.";
        return false;
    }

    int hdu_type(0);
    fits_movabs_hdu(fptr, hdu + 1, &hdu_type, &status);
    if (status || (hdu_type != BINARY_TBL && hdu_type != IMAGE_HDU)) {
        fits_close_file(fptr, &status);
        error = "HDU is not a binary table or image.";
        return false;
    }

    std::string pixtype;
    if (!ReadStringKey(fptr, "PIXTYPE", pixtype) || pixtype != "HEALPIX") {
        fits_close_file(fptr, &status);
        error = "HDU is not a HEALPix map.";
        return false;
    }

    if (!ReadLongKey(fptr, "NSIDE", info.nside) || info.nside <= 0) {
        fits_close_file(fptr, &status);
        error = "HEALPix map has no valid NSIDE.";
        return false;
    }

    std::string value;
    info.nested = ReadStringKey(fptr, "ORDERING", value) && value.rfind("NEST", 0) == 0;
    info.explicit_index = ReadStringKey(fptr, "INDXSCHM", value) && value == "EXPLICIT";
    if (!ReadStringKey(fptr, "COORDSYS", info.coordsys)) {
        info.coordsys = "";
    }
    if (!ReadStringKey(fptr, "OBJECT", info.object)) {
        info.object = info.explicit_index ? "PARTIAL" : "FULLSKY";
    }
    if (!ReadDoubleKey(fptr, "BAD_DATA", info.bad_data)) {
        info.bad_data = -1.6375e30;
    }
    if (!ReadDoubleKey(fptr, "BMAJ", info.bmaj_deg) || !ReadDoubleKey(fptr, "BMIN", info.bmin_deg)) {
        info.bmaj_deg = info.bmin_deg = std::nan("");
    }
    if (!ReadDoubleKey(fptr, "BPA", info.bpa_deg)) {
        info.bpa_deg = 0.0;
    }

    if (hdu_type == IMAGE_HDU) {
        // One map per row: the first axis is the HEALPix pixel index, the other axes enumerate the maps
        int naxis(0);
        long naxes[9] = {0};
        fits_get_img_dim(fptr, &naxis, &status);
        if (status || naxis < 1 || naxis > 9) {
            fits_close_file(fptr, &status);
            error = "HEALPix image HDU has no axes.";
            return false;
        }
        fits_get_img_size(fptr, naxis, naxes, &status);
        if (status || naxes[0] != info.npix() || info.explicit_index) {
            fits_close_file(fptr, &status);
            error = "HEALPix image HDU is not a full-sky map per row (NAXIS1 must be 12 NSIDE^2).";
            return false;
        }
        info.image_format = true;
        info.num_maps = 1;
        for (int i = 1; i < naxis; ++i) {
            info.num_maps *= naxes[i];
        }

        // WCS of the map axis: the first CTYPEn (n >= 2) present; files often number it as the cube's third axis
        info.map_axis_ctype.clear();
        for (int n = 2; n <= 4 && info.map_axis_ctype.empty(); ++n) {
            std::string ctype;
            char key[FLEN_KEYWORD];
            snprintf(key, sizeof(key), "CTYPE%d", n);
            if (ReadStringKey(fptr, key, ctype) && !ctype.empty()) {
                info.map_axis_ctype = ctype;
                snprintf(key, sizeof(key), "CUNIT%d", n);
                if (!ReadStringKey(fptr, key, info.map_axis_unit)) {
                    info.map_axis_unit = "";
                }
                snprintf(key, sizeof(key), "CRVAL%d", n);
                if (!ReadDoubleKey(fptr, key, info.map_axis_crval)) {
                    info.map_axis_crval = 0.0;
                }
                snprintf(key, sizeof(key), "CDELT%d", n);
                if (!ReadDoubleKey(fptr, key, info.map_axis_cdelt)) {
                    info.map_axis_cdelt = 1.0;
                }
                snprintf(key, sizeof(key), "CRPIX%d", n);
                if (!ReadDoubleKey(fptr, key, info.map_axis_crpix)) {
                    info.map_axis_crpix = 1.0;
                }
            }
        }
        if (!ReadStringKey(fptr, "SPECSYS", info.specsys)) {
            info.specsys = "";
        }
        if (!ReadDoubleKey(fptr, "RESTFRQ", info.rest_freq_hz) && !ReadDoubleKey(fptr, "RESTFREQ", info.rest_freq_hz)) {
            info.rest_freq_hz = std::nan("");
        }

        std::string bunit;
        if (!ReadStringKey(fptr, "BUNIT", bunit)) {
            bunit = "";
        }
        info.columns = {"MAP"};
        info.column_units = {bunit};
        fits_close_file(fptr, &status);
        return true;
    }

    // Column names and units
    int num_cols(0);
    fits_get_num_cols(fptr, &num_cols, &status);
    info.columns.clear();
    info.column_units.clear();
    for (int i = 1; i <= num_cols; ++i) {
        std::string name, unit;
        char key[FLEN_KEYWORD];
        snprintf(key, sizeof(key), "TTYPE%d", i);
        ReadStringKey(fptr, key, name);
        snprintf(key, sizeof(key), "TUNIT%d", i);
        if (!ReadStringKey(fptr, key, unit)) {
            unit = "";
        }
        if (info.explicit_index && (name == "PIXEL" || name == "PIX")) {
            continue; // index column, not data
        }
        info.columns.push_back(name);
        info.column_units.push_back(unit);
    }

    fits_close_file(fptr, &status);
    return !info.columns.empty();
}

// ---------------------------------------------------------------------------
// Construction

CartaHealpixImage::CartaHealpixImage(const std::string& filename, unsigned int hdu, const HealpixProjectionOptions& options)
    : casacore::ImageInterface<float>(),
      _filename(filename),
      _hdu(hdu),
      _options(options),
      _grid_frame(casacore::MDirection::ICRS),
      _map_frame(casacore::MDirection::ICRS),
      _needs_rotation(false) {
    std::string error;
    if (!ReadMapInfo(filename, hdu, _info, error)) {
        throw(casacore::AipsError(error));
    }

    _healpix = std::make_shared<T_Healpix_Base<int64_t>>(_info.nside, _info.nested ? NEST : RING, SET_NSIDE);
    _map_frame = _info.frame();
    if (_info.coordsys.empty()) {
        spdlog::warn("HEALPix map {} has no COORDSYS header; assuming celestial (ICRS).", filename);
    }

    ReadMap();
    SetUpCoordinateSystem();
    SetUpFrameRotation();
}

CartaHealpixImage::CartaHealpixImage(const CartaHealpixImage& other)
    : casacore::ImageInterface<float>(other),
      _filename(other._filename),
      _hdu(other._hdu),
      _options(other._options),
      _info(other._info),
      _column_names(other._column_names),
      _source(other._source),
      _healpix(other._healpix),
      _shape(other._shape),
      _direction_coord(other._direction_coord),
      _grid_frame(other._grid_frame),
      _map_frame(other._map_frame),
      _needs_rotation(other._needs_rotation),
      _rotation(other._rotation) {}

void CartaHealpixImage::ReadMap() {
    if (_info.image_format) {
        ReadImageFormatMaps();
        return;
    }

    // Select columns
    std::vector<int> selected;
    if (_options.columns.empty()) {
        for (size_t i = 0; i < _info.columns.size(); ++i) {
            selected.push_back(i);
        }
    } else {
        for (auto& requested : _options.columns) {
            auto name = UpperTrim(requested);
            bool found(false);
            for (size_t i = 0; i < _info.columns.size(); ++i) {
                if (UpperTrim(_info.columns[i]) == name) {
                    selected.push_back(i);
                    found = true;
                    break;
                }
            }
            if (!found) {
                throw(casacore::AipsError("HEALPix column " + requested + " not found."));
            }
        }
    }

    fitsfile* fptr(nullptr);
    int status(0);
    fits_open_file(&fptr, _filename.c_str(), READONLY, &status);
    int hdu_type(0);
    fits_movabs_hdu(fptr, _hdu + 1, &hdu_type, &status);
    if (status) {
        throw(casacore::AipsError("Cannot open HEALPix HDU."));
    }

    long num_rows(0);
    fits_get_num_rows(fptr, &num_rows, &status);
    const long npix = _info.npix();

    // Explicit (partial-sky) pixel indices
    std::vector<int64_t> pixel_index;
    long num_values(0);
    if (_info.explicit_index) {
        int pix_col(0);
        char pix_name[] = "PIXEL";
        fits_get_colnum(fptr, CASEINSEN, pix_name, &pix_col, &status);
        if (status) {
            status = 0;
            char pix_name2[] = "PIX";
            fits_get_colnum(fptr, CASEINSEN, pix_name2, &pix_col, &status);
        }
        if (status) {
            fits_close_file(fptr, &status);
            throw(casacore::AipsError("HEALPix partial map has no PIXEL column."));
        }
        int typecode(0);
        long repeat(1), width(0);
        fits_get_coltype(fptr, pix_col, &typecode, &repeat, &width, &status);
        num_values = num_rows * repeat;
        pixel_index.resize(num_values);
        fits_read_col(fptr, TLONGLONG, pix_col, 1, 1, num_values, nullptr, pixel_index.data(), nullptr, &status);
        if (status) {
            fits_close_file(fptr, &status);
            throw(casacore::AipsError("Error reading HEALPix PIXEL column."));
        }
    }

    _column_names.clear();
    std::vector<std::shared_ptr<std::vector<float>>> maps;
    for (int idx : selected) {
        const std::string& name = _info.columns[idx];
        int colnum(0);
        std::vector<char> cname(name.begin(), name.end());
        cname.push_back('\0');
        fits_get_colnum(fptr, CASEINSEN, cname.data(), &colnum, &status);
        if (status) {
            fits_close_file(fptr, &status);
            throw(casacore::AipsError("HEALPix column " + name + " not found in table."));
        }

        int typecode(0);
        long repeat(1), width(0);
        fits_get_coltype(fptr, colnum, &typecode, &repeat, &width, &status);
        long column_values = num_rows * repeat;

        auto map = std::make_shared<std::vector<float>>(npix, FLOAT_NAN);
        if (_info.explicit_index) {
            if (column_values != num_values) {
                fits_close_file(fptr, &status);
                throw(casacore::AipsError("HEALPix PIXEL and data columns have different lengths."));
            }
            std::vector<float> values(column_values);
            fits_read_col(fptr, TFLOAT, colnum, 1, 1, column_values, nullptr, values.data(), nullptr, &status);
            for (long i = 0; i < column_values; ++i) {
                auto pix = pixel_index[i];
                if (pix >= 0 && pix < npix) {
                    (*map)[pix] = values[i];
                }
            }
        } else {
            long nread = std::min(column_values, npix);
            if (column_values != npix) {
                spdlog::warn("HEALPix column {} has {} values, expected {} for NSIDE={}.", name, column_values, npix, _info.nside);
            }
            fits_read_col(fptr, TFLOAT, colnum, 1, 1, nread, nullptr, map->data(), nullptr, &status);
        }
        if (status) {
            fits_close_file(fptr, &status);
            throw(casacore::AipsError("Error reading HEALPix column " + name + "."));
        }

        // Replace UNSEEN / BAD_DATA with NaN
        const float bad = static_cast<float>(_info.bad_data);
        for (auto& v : *map) {
            if (v == bad || v < -1.0e30f) {
                v = FLOAT_NAN;
            }
        }

        _column_names.push_back(name);
        maps.push_back(map);
    }

    fits_close_file(fptr, &status);
    _source = std::make_shared<InMemoryMapSource>(std::move(maps));
}

void CartaHealpixImage::ReadImageFormatMaps() {
    // Maps are rows of the image HDU; they are read lazily. An uncompressed float file is memory-mapped so that
    // a cube of hundreds of maps costs no memory of its own; anything else goes through cfitsio with a small cache.
    _column_names = _info.columns;
    const float bad = static_cast<float>(_info.bad_data);

    fitsfile* fptr(nullptr);
    int status(0);
    fits_open_file(&fptr, _filename.c_str(), READONLY, &status);
    int hdu_type(0);
    fits_movabs_hdu(fptr, _hdu + 1, &hdu_type, &status);
    if (status) {
        throw(casacore::AipsError("Cannot open HEALPix image HDU."));
    }
    int bitpix(0);
    fits_get_img_type(fptr, &bitpix, &status);
    double bscale(1.0), bzero(0.0);
    if (!ReadDoubleKey(fptr, "BSCALE", bscale)) {
        bscale = 1.0;
    }
    if (!ReadDoubleKey(fptr, "BZERO", bzero)) {
        bzero = 0.0;
    }
    LONGLONG headstart(0), datastart(0), dataend(0);
    ffghadll(fptr, &headstart, &datastart, &dataend, &status); // byte offsets of this HDU within the file
    int compressed = fits_is_compressed_image(fptr, &status);
    fits_close_file(fptr, &status);

    bool can_map = (status == 0) && !compressed && !IsCompressedFilename(_filename) && (bitpix == FLOAT_IMG || bitpix == DOUBLE_IMG);
    if (can_map) {
        try {
            _source =
                std::make_shared<MappedImageMapSource>(_filename, datastart, bitpix, _info.npix(), _info.num_maps, bscale, bzero, bad);
            return;
        } catch (const casacore::AipsError& err) {
            spdlog::debug("HEALPix image {}: {} Falling back to cfitsio reads.", _filename, err.getMesg());
        }
    }
    _source = std::make_shared<CfitsioRowMapSource>(_filename, _hdu, _info.npix(), _info.num_maps, bad);
}

size_t CartaHealpixImage::NumMaps() const {
    return _source ? _source->NumMaps() : 0;
}

void CartaHealpixImage::SetUpCoordinateSystem() {
    casacore::CoordinateSystem coord_sys;
    int width(0), height(0);

    if (_options.match_coord_sys && _options.match_coord_sys->hasDirectionCoordinate()) {
        // Reproject onto another image's grid: keep only its direction coordinate and pixel axes
        casacore::CoordinateSystem matched(*_options.match_coord_sys);
        auto dir_axes = matched.directionAxesNumbers();
        int x_axis = std::min(dir_axes(0), dir_axes(1));
        int y_axis = std::max(dir_axes(0), dir_axes(1));
        if (_options.match_shape.size() <= static_cast<size_t>(y_axis)) {
            throw(casacore::AipsError("Matched image shape does not cover its direction axes."));
        }
        width = _options.match_shape(x_axis);
        height = _options.match_shape(y_axis);

        // Remove every world axis that is not a direction axis (removes its pixel axis too).
        // The replacement value must be valid for the axis (e.g. a Stokes type), so use its reference value.
        int dir_coord = matched.directionCoordinateNumber();
        for (int axis = matched.nWorldAxes() - 1; axis >= 0; --axis) {
            int coord(-1), axis_in_coord(-1);
            matched.findWorldAxis(coord, axis_in_coord, axis);
            if (coord != dir_coord) {
                double replacement = matched.referenceValue()(axis);
                if (!matched.removeWorldAxis(axis, replacement)) {
                    throw(casacore::AipsError("Cannot remove axis from the matched image coordinate system: " + matched.errorMessage()));
                }
            }
        }
        if (matched.nPixelAxes() != 2 || matched.nWorldAxes() != 2) {
            throw(casacore::AipsError("Cannot extract a 2D direction grid from the matched image."));
        }
        if (matched.directionAxesNumbers()(0) > matched.directionAxesNumbers()(1)) {
            // Latitude is the first pixel axis in the matched image; our x axis follows that order too
            casacore::Vector<casacore::Int> world_order(2), pixel_order(2);
            world_order(0) = 0;
            world_order(1) = 1;
            pixel_order(0) = 1;
            pixel_order(1) = 0;
            matched.transpose(world_order, pixel_order);
        }
        int new_width = _options.width > 0 ? _options.width : width;
        int new_height = _options.height > 0 ? _options.height : height;
        if (_options.extend_sky && _options.width <= 0 && _options.height <= 0) {
            // Show the map around the image: enlarge the grid, centred on the image, out to where the projection stops
            // being valid (the hemisphere of a SIN image, the whole sphere of a ZEA image) at the image's pixel scale
            auto extent = ExtendedGridSize(matched.directionCoordinate(), width, height);
            new_width = extent.first;
            new_height = extent.second;
        }
        if (new_width != width || new_height != height) {
            // A different grid size keeps the image's WCS and centres its grid in ours, so overlays still align
            casacore::DirectionCoordinate dir = matched.directionCoordinate();
            casacore::Vector<casacore::Double> refpix = dir.referencePixel();
            refpix(0) += (new_width - width) / 2.0;
            refpix(1) += (new_height - height) / 2.0;
            dir.setReferencePixel(refpix);
            matched.replaceCoordinate(dir, matched.directionCoordinateNumber());
            width = new_width;
            height = new_height;
        }
        coord_sys = matched;
        _direction_coord = matched.directionCoordinate();
        _grid_frame = _direction_coord.directionType();
    } else {
        // Build our own projection
        casacore::Projection::Type proj_type = casacore::Projection::type(casacore::String(UpperTrim(_options.projection)));
        if (proj_type == casacore::Projection::N_PROJ) {
            throw(casacore::AipsError("Unknown projection " + _options.projection + " for HEALPix map."));
        }

        bool frame_ok(false);
        _grid_frame = _options.frame.empty() ? _map_frame : FrameFromCode(_options.frame, frame_ok);
        if (!_options.frame.empty() && !frame_ok) {
            throw(casacore::AipsError("Unknown coordinate frame " + _options.frame + " for HEALPix map."));
        }

        double scale =
            std::isnan(_options.pixel_scale_deg) || _options.pixel_scale_deg <= 0 ? _info.resolution_deg() : _options.pixel_scale_deg;
        double lon0 = std::isnan(_options.center_lon_deg) ? 0.0 : _options.center_lon_deg;
        double lat0 = std::isnan(_options.center_lat_deg) ? 0.0 : _options.center_lat_deg;

        casacore::Matrix<double> xform(2, 2, 0.0);
        xform.diagonal() = 1.0;

        // Grid size: as far as the projection reaches from its reference point along each axis (the whole ellipse of
        // a Mollweide or Aitoff map, the full hemisphere of a SIN map, the whole sphere of a ZEA map), plus one pixel
        // of margin, with the reference pixel at the centre; never more than 360 degrees across
        width = _options.width;
        height = _options.height;
        if (width <= 0 || height <= 0) {
            const int max_half = std::min(MAX_GRID_SIDE / 2, static_cast<int>(std::ceil(180.0 / scale)));
            casacore::DirectionCoordinate probe(_grid_frame, casacore::Projection(proj_type), lon0 * DEG2RAD, lat0 * DEG2RAD,
                -scale * DEG2RAD, scale * DEG2RAD, xform, 0.0, 0.0);
            auto reach = ValidReach(probe, 0.0, 0.0, max_half);
            if (width <= 0) {
                width = 2 * (reach.first + 1) + 1;
            }
            if (height <= 0) {
                height = 2 * (reach.second + 1) + 1;
            }
        }
        // Reference pixel at the image centre (0-based), longitude increasing to the left
        _direction_coord = casacore::DirectionCoordinate(_grid_frame, casacore::Projection(proj_type), lon0 * DEG2RAD, lat0 * DEG2RAD,
            -scale * DEG2RAD, scale * DEG2RAD, xform, (width - 1) / 2.0, (height - 1) / 2.0);
        coord_sys.addCoordinate(_direction_coord);
    }

    // Third axis: the maps
    const size_t ncol = NumMaps();
    if (_info.image_format && ncol > 1) {
        if (!AddMapAxisCoordinate(coord_sys)) {
            casacore::Vector<casacore::String> names(1, "MAP"), units(1, "");
            casacore::Vector<casacore::Double> refval(1, 0.0), inc(1, 1.0), refpix(1, 0.0);
            casacore::Matrix<casacore::Double> pc(1, 1, 1.0);
            coord_sys.addCoordinate(casacore::LinearCoordinate(names, units, refval, inc, pc, refpix));
        }
        _shape = casacore::IPosition(3, width, height, ncol);
    } else if (ncol > 1) {
        casacore::Vector<casacore::Int> stokes(ncol);
        bool all_stokes(true);
        for (size_t i = 0; i < ncol; ++i) {
            auto type = ColumnStokes(_column_names[i]);
            stokes(i) = type;
            all_stokes &= (type != casacore::Stokes::Undefined);
        }
        // A Stokes coordinate must not repeat a type
        for (size_t i = 0; all_stokes && i < ncol; ++i) {
            for (size_t j = i + 1; j < ncol; ++j) {
                all_stokes &= (stokes(i) != stokes(j));
            }
        }

        if (all_stokes) {
            coord_sys.addCoordinate(casacore::StokesCoordinate(stokes));
        } else {
            casacore::Vector<casacore::String> names(1, "COLUMN"), units(1, "");
            casacore::Vector<casacore::Double> refval(1, 0.0), inc(1, 1.0), refpix(1, 0.0);
            casacore::Matrix<casacore::Double> pc(1, 1, 1.0);
            coord_sys.addCoordinate(casacore::LinearCoordinate(names, units, refval, inc, pc, refpix));
        }
        _shape = casacore::IPosition(3, width, height, ncol);
    } else {
        _shape = casacore::IPosition(2, width, height);
    }

    setCoordinateInfo(coord_sys);

    // Brightness unit from the first column, if casacore understands it
    if (!_info.column_units.empty()) {
        casacore::UnitMap::addFITS();
        std::string unit = _info.column_units[0];
        if (!unit.empty() && casacore::UnitVal::check(unit)) {
            setUnits(casacore::Unit(unit));
        }
    }

    // Restoring beam from BMAJ/BMIN/BPA (degrees)
    if (std::isfinite(_info.bmaj_deg) && std::isfinite(_info.bmin_deg) && _info.bmaj_deg > 0 && _info.bmin_deg > 0) {
        casacore::ImageInfo info = imageInfo();
        info.setRestoringBeam(casacore::GaussianBeam(casacore::Quantity(_info.bmaj_deg, "deg"), casacore::Quantity(_info.bmin_deg, "deg"),
            casacore::Quantity(std::isfinite(_info.bpa_deg) ? _info.bpa_deg : 0.0, "deg")));
        setImageInfo(info);
    }
}

bool CartaHealpixImage::AddMapAxisCoordinate(casacore::CoordinateSystem& coord_sys) {
    // Spectral axis for FREQ-like map axes, a linear axis for any other named axis
    if (!_info.has_map_axis()) {
        return false;
    }
    std::string ctype = UpperTrim(_info.map_axis_ctype);
    std::string base = ctype.substr(0, ctype.find('-'));
    casacore::UnitMap::addFITS();
    const double refpix = _info.map_axis_crpix - 1.0; // 0-based

    if (base == "FREQ") {
        double factor(1.0);
        if (!_info.map_axis_unit.empty() && casacore::UnitVal::check(_info.map_axis_unit)) {
            factor = casacore::Quantity(1.0, _info.map_axis_unit).getValue("Hz");
        }
        casacore::MFrequency::Types freq_type(casacore::MFrequency::TOPO);
        if (!_info.specsys.empty()) {
            casacore::MFrequency::Types parsed;
            if (casacore::MFrequency::getType(parsed, _info.specsys)) {
                freq_type = parsed;
            }
        }
        double rest_freq = std::isfinite(_info.rest_freq_hz) ? _info.rest_freq_hz : 0.0;
        coord_sys.addCoordinate(
            casacore::SpectralCoordinate(freq_type, _info.map_axis_crval * factor, _info.map_axis_cdelt * factor, refpix, rest_freq));
        return true;
    }

    casacore::Vector<casacore::String> names(1, base.empty() ? "MAP" : base), units(1, "");
    if (!_info.map_axis_unit.empty() && casacore::UnitVal::check(_info.map_axis_unit)) {
        units(0) = _info.map_axis_unit;
    }
    casacore::Vector<casacore::Double> refval(1, _info.map_axis_crval), inc(1, _info.map_axis_cdelt), refpix_v(1, refpix);
    casacore::Matrix<casacore::Double> pc(1, 1, 1.0);
    coord_sys.addCoordinate(casacore::LinearCoordinate(names, units, refval, inc, pc, refpix_v));
    return true;
}

std::pair<int, int> CartaHealpixImage::ValidReach(const casacore::DirectionCoordinate& direction, double cx, double cy, int max_half) {
    // Largest number of pixels from (cx, cy) along each grid axis at which the projection is still valid: the pixel
    // converts to a sky position, and that position converts back to the same pixel (a plate carree grid, say,
    // converts pixels beyond 180 degrees of longitude to positions that belong to pixels on its other side).
    // Found by bisection, so at most max_half.
    auto on_sky = [&](double x, double y) {
        casacore::Vector<casacore::Double> pixel(2), world(2), back(2);
        pixel(0) = x;
        pixel(1) = y;
        if (!direction.toWorld(world, pixel) || !std::isfinite(world(0)) || !std::isfinite(world(1))) {
            return false;
        }
        if (!direction.toPixel(back, world) || !std::isfinite(back(0)) || !std::isfinite(back(1))) {
            return false;
        }
        return std::abs(back(0) - x) < 0.5 && std::abs(back(1) - y) < 0.5;
    };
    auto reach = [&](double sx, double sy) { // along the direction (sx, sy) from the centre
        int lo(0), hi(max_half);
        while (lo < hi) {
            int mid = (lo + hi + 1) / 2;
            if (on_sky(cx + sx * mid, cy + sy * mid)) {
                lo = mid;
            } else {
                hi = mid - 1;
            }
        }
        return lo;
    };
    return {std::max(reach(1, 0), reach(-1, 0)), std::max(reach(0, 1), reach(0, -1))};
}

std::pair<int, int> CartaHealpixImage::ExtendedGridSize(const casacore::DirectionCoordinate& direction, int width, int height) {
    // The grid, centred on the image, enlarged to everywhere the image's projection is valid (the hemisphere of a
    // SIN image, the whole sphere of a ZEA image, as far as the cap allows for a TAN image): all the sky the map
    // can show around the image, at the image's pixel scale
    const double cx = (width - 1) / 2.0, cy = (height - 1) / 2.0;
    auto reach = ValidReach(direction, cx, cy, MAX_GRID_SIDE / 2);
    int new_width = std::max(width, std::min(MAX_GRID_SIDE, 2 * (reach.first + 1) + 1));
    int new_height = std::max(height, std::min(MAX_GRID_SIDE, 2 * (reach.second + 1) + 1));
    // Keep the image's pixels on the new grid's integer positions: the offset (new - old) / 2 must be a whole number
    if ((new_width - width) % 2) {
        new_width = (new_width < MAX_GRID_SIDE) ? new_width + 1 : new_width - 1;
    }
    if ((new_height - height) % 2) {
        new_height = (new_height < MAX_GRID_SIDE) ? new_height + 1 : new_height - 1;
    }
    return {new_width, new_height};
}

void CartaHealpixImage::SetUpFrameRotation() {
    _needs_rotation = (_grid_frame != _map_frame);
    if (!_needs_rotation) {
        return;
    }

    // Build the 3x3 rotation from the grid frame to the map frame by converting the basis vectors.
    // This is exact for frames related by a pure rotation (ICRS/J2000, Galactic, ecliptic, supergalactic).
    casacore::MDirection::Ref grid_ref(_grid_frame), map_ref(_map_frame);
    casacore::MDirection::Convert converter(grid_ref, map_ref);
    _rotation.resize(3, 3);
    for (int i = 0; i < 3; ++i) {
        casacore::Vector<casacore::Double> basis(3, 0.0);
        basis(i) = 1.0;
        casacore::MDirection in(casacore::MVDirection(basis), _grid_frame);
        casacore::Vector<casacore::Double> out = converter(in).getValue().getValue();
        for (int j = 0; j < 3; ++j) {
            _rotation(j, i) = out(j);
        }
    }

    // Sanity check on an arbitrary direction; warn if the frames are not related by a rotation
    casacore::Vector<casacore::Double> test(3);
    test(0) = 0.3;
    test(1) = -0.5;
    test(2) = 0.812403840463596; // unit vector
    casacore::MDirection in(casacore::MVDirection(test), _grid_frame);
    casacore::Vector<casacore::Double> expected = converter(in).getValue().getValue();
    double err(0.0);
    for (int j = 0; j < 3; ++j) {
        double rotated = _rotation(j, 0) * test(0) + _rotation(j, 1) * test(1) + _rotation(j, 2) * test(2);
        err = std::max(err, std::abs(rotated - expected(j)));
    }
    if (err > 1.0e-6) {
        spdlog::warn("HEALPix frame conversion {} -> {} is not a pure rotation (error {:.2e} rad); positions may be slightly off.",
            std::string(casacore::MDirection::showType(_grid_frame)), std::string(casacore::MDirection::showType(_map_frame)), err);
    }
}

// ---------------------------------------------------------------------------
// casacore::ImageInterface

casacore::String CartaHealpixImage::imageType() const {
    return "CartaHealpixImage";
}

casacore::String CartaHealpixImage::name(bool stripPath) const {
    if (stripPath) {
        auto pos = _filename.find_last_of('/');
        return pos == std::string::npos ? _filename : _filename.substr(pos + 1);
    }
    return _filename;
}

casacore::IPosition CartaHealpixImage::shape() const {
    return _shape;
}

casacore::Bool CartaHealpixImage::ok() const {
    return _source && _source->NumMaps() > 0;
}

casacore::ImageInterface<float>* CartaHealpixImage::cloneII() const {
    return new CartaHealpixImage(*this);
}

void CartaHealpixImage::doPutSlice(const casacore::Array<float>&, const casacore::IPosition&, const casacore::IPosition&) {
    throw(casacore::AipsError("CartaHealpixImage is read-only."));
}

const casacore::LatticeRegion* CartaHealpixImage::getRegionPtr() const {
    return nullptr;
}

void CartaHealpixImage::resize(const casacore::TiledShape&) {
    throw(casacore::AipsError("CartaHealpixImage cannot be resized."));
}

casacore::uInt CartaHealpixImage::advisedMaxPixels() const {
    return CURSOR_SIZE * CURSOR_SIZE;
}

casacore::IPosition CartaHealpixImage::doNiceCursorShape(casacore::uInt) const {
    casacore::IPosition cursor(_shape.size(), 1);
    cursor(0) = std::min<ssize_t>(_shape(0), CURSOR_SIZE);
    cursor(1) = std::min<ssize_t>(_shape(1), CURSOR_SIZE);
    return cursor;
}

casacore::Bool CartaHealpixImage::isMasked() const {
    return false; // blanks are NaN
}

casacore::Bool CartaHealpixImage::hasPixelMask() const {
    return false;
}

casacore::Bool CartaHealpixImage::doGetMaskSlice(casacore::Array<bool>& buffer, const casacore::Slicer& section) {
    buffer.resize(section.length());
    buffer = true;
    return false;
}

casacore::Bool CartaHealpixImage::doGetSlice(casacore::Array<float>& buffer, const casacore::Slicer& section) {
    const casacore::IPosition blc = section.start();
    const casacore::IPosition inc = section.stride();
    const casacore::IPosition length = section.length();
    const size_t nx = length(0);
    const size_t ny = length(1);
    const size_t ncol = length.size() > 2 ? length(2) : 1;

    if (!buffer.shape().isEqual(length)) {
        buffer.resize(length);
    }
    if (nx * ny == 0) {
        return false;
    }

    bool delete_it;
    float* data = buffer.getStorage(delete_it);
    Sample(blc(0), inc(0), nx, blc(1), inc(1), ny, length.size() > 2 ? blc(2) : 0, length.size() > 2 ? inc(2) : 1, ncol, data);
    buffer.putStorage(data, delete_it);
    // As for CartaFitsImage, the return value is used by FileLoader::GetSlice as a success flag
    return true;
}

void CartaHealpixImage::Sample(double x0, double dx, size_t nx, double y0, double dy, size_t ny, size_t column0, size_t column_stride,
    size_t ncol, float* data) const {
    const size_t npoints = nx * ny;
    if (npoints == 0) {
        return;
    }

    // The points are processed in chunks, in parallel, so that a full-resolution plane of a large grid (tens of
    // millions of pixels when matched to a survey image) never needs its coordinates held in memory at once.
    const size_t chunk_size(1 << 16);
    const size_t num_chunks = (npoints + chunk_size - 1) / chunk_size;

#pragma omp parallel for schedule(dynamic)
    for (size_t chunk = 0; chunk < num_chunks; ++chunk) {
        const size_t chunk_start = chunk * chunk_size;
        const size_t nchunk = std::min(chunk_size, npoints - chunk_start);

        // Grid pixel coordinates (0-based) of the chunk's points, x fastest
        std::vector<double> x(nchunk), y(nchunk);
        for (size_t k = 0; k < nchunk; ++k) {
            const size_t p = chunk_start + k;
            x[k] = x0 + (p % nx) * dx;
            y[k] = y0 + (p / nx) * dy;
        }

        HealpixLookupPlan plan;
        PlanLookup(x.data(), y.data(), nchunk, plan);
        for (size_t c = 0; c < ncol; ++c) {
            LookupPlanned(plan, column0 + c * column_stride, data + c * npoints + chunk_start);
        }
    }
}

size_t CartaHealpixImage::Npix() const {
    return _info.npix();
}

void CartaHealpixImage::ReadMapRange(size_t map, int64_t first, size_t count, float* out) const {
    if (!_source || map >= _source->NumMaps()) {
        std::fill(out, out + count, FLOAT_NAN);
        return;
    }
    _source->ReadRange(map, first, count, out);
}

void CartaHealpixImage::ToMapDirections(
    const double* x, const double* y, size_t npoints, std::vector<double>& lon, std::vector<double>& lat, std::vector<bool>& valid) const {
    // The direction coordinate's own pixel-axis order (the image may have latitude as its first axis when matched)
    const auto& coord_sys = coordinates();
    const int dir_coord = coord_sys.directionCoordinateNumber();
    const casacore::DirectionCoordinate& direction = coord_sys.directionCoordinate();
    auto dir_pixel_axes = coord_sys.pixelAxes(dir_coord); // image axis of the coordinate's pixel axes 0 (lon) and 1 (lat)
    const int lon_row = (dir_pixel_axes(0) == 0) ? 0 : 1;
    const int lat_row = 1 - lon_row;

    auto units = direction.worldAxisUnits();
    const double lon_factor = casacore::Quantity(1.0, units(0)).getValue("rad");
    const double lat_factor = casacore::Quantity(1.0, units(1)).getValue("rad");
    const double max_lat = PI / 2.0 + 1.0e-9;

    // Each call converts with its own copy of the direction coordinate (casacore's is not thread-safe)
    casacore::DirectionCoordinate direction_copy(direction);
    casacore::Matrix<casacore::Double> pixel(2, npoints, 0.0);
    for (size_t k = 0; k < npoints; ++k) {
        pixel(lon_row, k) = x[k];
        pixel(lat_row, k) = y[k];
    }
    lon.assign(npoints, 0.0);
    lat.assign(npoints, 0.0);
    valid.assign(npoints, false);

    // Convert to world (radians), in blocks: casacore's DirectionCoordinate::toWorldMany returns garbage for the
    // points after the first off-sky pixel of a batch, so any block with a failure is redone point by point
    const size_t block_size(4096);
    casacore::Vector<casacore::Double> one_pixel(2), one_world(2);
    for (size_t block_start = 0; block_start < npoints; block_start += block_size) {
        const size_t n = std::min(npoints, block_start + block_size) - block_start;
        casacore::Matrix<casacore::Double> block_pixel = pixel(casacore::Slice(0, 2), casacore::Slice(block_start, n));
        casacore::Matrix<casacore::Double> world;
        casacore::Vector<casacore::Bool> failures;

        if (direction_copy.toWorldMany(world, block_pixel, failures)) {
            for (size_t i = 0; i < n; ++i) {
                lon[block_start + i] = world(0, i) * lon_factor;
                lat[block_start + i] = world(1, i) * lat_factor;
                valid[block_start + i] = true;
            }
        } else {
            for (size_t i = 0; i < n; ++i) {
                one_pixel(0) = block_pixel(0, i);
                one_pixel(1) = block_pixel(1, i);
                if (direction_copy.toWorld(one_world, one_pixel)) {
                    lon[block_start + i] = one_world(0) * lon_factor;
                    lat[block_start + i] = one_world(1) * lat_factor;
                    valid[block_start + i] = true;
                }
            }
        }
    }

    for (size_t k = 0; k < npoints; ++k) {
        if (!valid[k]) {
            continue;
        }
        if (!std::isfinite(lon[k]) || !std::isfinite(lat[k]) || std::abs(lat[k]) > max_lat) {
            valid[k] = false;
        } else {
            lat[k] = std::max(-PI / 2.0, std::min(PI / 2.0, lat[k]));
        }
    }

    if (_needs_rotation) {
        for (size_t k = 0; k < npoints; ++k) {
            if (!valid[k]) {
                continue;
            }
            double cl = std::cos(lat[k]);
            double v0 = cl * std::cos(lon[k]), v1 = cl * std::sin(lon[k]), v2 = std::sin(lat[k]);
            double r0 = _rotation(0, 0) * v0 + _rotation(0, 1) * v1 + _rotation(0, 2) * v2;
            double r1 = _rotation(1, 0) * v0 + _rotation(1, 1) * v1 + _rotation(1, 2) * v2;
            double r2 = _rotation(2, 0) * v0 + _rotation(2, 1) * v1 + _rotation(2, 2) * v2;
            lon[k] = std::atan2(r1, r0);
            lat[k] = std::asin(std::max(-1.0, std::min(1.0, r2)));
        }
    }
}

void CartaHealpixImage::PlanLookup(const double* x, const double* y, size_t npoints, HealpixLookupPlan& plan) const {
    std::vector<double> lon, lat;
    std::vector<bool> valid;
    ToMapDirections(x, y, npoints, lon, lat, valid);

    plan.npoints = npoints;
    plan.bilinear = _options.bilinear;
    const size_t neighbours = plan.neighbours();
    plan.point_pixels.assign(npoints * neighbours, -1);
    plan.point_weights.clear();
    if (plan.bilinear) {
        plan.point_weights.assign(npoints * 4, 0.0f);
    }
    plan.pixels.clear();

    // HEALPix pixel (or 4 neighbours and weights) of every point on the sky
    std::vector<int64_t> wanted;
    wanted.reserve(npoints * neighbours);
    for (size_t p = 0; p < npoints; ++p) {
        if (!valid[p]) {
            continue;
        }
        double theta = PI / 2.0 - lat[p];
        double phi = lon[p];
        if (phi < 0.0) {
            phi += 2.0 * PI;
        }
        pointing ptg(theta, phi);
        if (plan.bilinear) {
            fix_arr<int64_t, 4> pix;
            fix_arr<double, 4> wgt;
            _healpix->get_interpol(ptg, pix, wgt);
            for (int n = 0; n < 4; ++n) {
                wanted.push_back(pix[n]);
                plan.point_weights[p * 4 + n] = static_cast<float>(wgt[n]);
            }
        } else {
            wanted.push_back(_healpix->ang2pix(ptg));
        }
    }

    // Distinct pixels, ascending, and each point's indices into them
    plan.pixels = wanted;
    std::sort(plan.pixels.begin(), plan.pixels.end());
    plan.pixels.erase(std::unique(plan.pixels.begin(), plan.pixels.end()), plan.pixels.end());
    size_t next(0);
    for (size_t p = 0; p < npoints; ++p) {
        if (!valid[p]) {
            continue;
        }
        for (size_t n = 0; n < neighbours; ++n) {
            auto found = std::lower_bound(plan.pixels.begin(), plan.pixels.end(), wanted[next++]);
            plan.point_pixels[p * neighbours + n] = static_cast<int32_t>(found - plan.pixels.begin());
        }
    }
}

void CartaHealpixImage::LookupPlanned(const HealpixLookupPlan& plan, size_t map, float* out) const {
    if (!_source || map >= _source->NumMaps() || plan.pixels.empty()) {
        std::fill(out, out + plan.npoints, FLOAT_NAN);
        return;
    }
    std::vector<float> values(plan.pixels.size());
    _source->ReadPixels(map, plan.pixels.data(), plan.pixels.size(), values.data());

    if (plan.bilinear) {
        for (size_t p = 0; p < plan.npoints; ++p) {
            if (plan.point_pixels[p * 4] < 0) {
                out[p] = FLOAT_NAN;
                continue;
            }
            double sum(0.0), wsum(0.0);
            for (int n = 0; n < 4; ++n) {
                float v = values[plan.point_pixels[p * 4 + n]];
                if (std::isfinite(v)) {
                    sum += plan.point_weights[p * 4 + n] * v;
                    wsum += plan.point_weights[p * 4 + n];
                }
            }
            out[p] = wsum > 0.0 ? static_cast<float>(sum / wsum) : FLOAT_NAN;
        }
    } else {
        for (size_t p = 0; p < plan.npoints; ++p) {
            int32_t index = plan.point_pixels[p];
            out[p] = index < 0 ? FLOAT_NAN : values[index];
        }
    }
}
