/* This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-backend
   Copyright 2018- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
   Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
   SPDX-License-Identifier: GPL-3.0-or-later
*/

//# CartaHealpixImage.cc : HEALPix map presented as a casacore::ImageInterface via on-the-fly reprojection

#include "CartaHealpixImage.h"

#include <algorithm>
#include <cmath>
#include <limits>

#include <fitsio.h>

#include <casacore/casa/Arrays/Slice.h>
#include <casacore/casa/Exceptions/Error.h>
#include <casacore/casa/Quanta/MVDirection.h>
#include <casacore/casa/Quanta/Quantum.h>
#include <casacore/casa/Quanta/UnitMap.h>
#include <casacore/casa/Quanta/UnitVal.h>
#include <casacore/coordinates/Coordinates/LinearCoordinate.h>
#include <casacore/coordinates/Coordinates/Projection.h>
#include <casacore/coordinates/Coordinates/StokesCoordinate.h>
#include <casacore/measures/Measures/MCDirection.h>
#include <casacore/measures/Measures/MeasConvert.h>
#include <casacore/measures/Measures/Stokes.h>

#include <healpix_base.h>
#include <pointing.h>

#include "Logger/Logger.h"
#include "Util/Nan.h"

using namespace carta;

namespace {

const double PI = M_PI;
const double DEG2RAD = M_PI / 180.0;
const unsigned int CURSOR_SIZE = 512;

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

bool IsZenithal(casacore::Projection::Type type) {
    switch (type) {
        case casacore::Projection::AZP:
        case casacore::Projection::SZP:
        case casacore::Projection::TAN:
        case casacore::Projection::SIN:
        case casacore::Projection::STG:
        case casacore::Projection::ARC:
        case casacore::Projection::ZPN:
        case casacore::Projection::ZEA:
        case casacore::Projection::AIR:
            return true;
        default:
            return false;
    }
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

} // namespace

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
    // HEALPix maps are binary tables in an extension; check the first few extensions.
    fitsfile* fptr(nullptr);
    int status(0);
    fits_open_file(&fptr, filename.c_str(), READONLY, &status);
    if (status) {
        return false;
    }
    int num_hdus(0);
    fits_get_num_hdus(fptr, &num_hdus, &status);
    fits_close_file(fptr, &status);

    for (int hdu = 1; hdu < std::min(num_hdus, 4); ++hdu) {
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
    if (status || hdu_type != BINARY_TBL) {
        fits_close_file(fptr, &status);
        error = "HDU is not a binary table.";
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
      _maps(other._maps),
      _healpix(other._healpix),
      _shape(other._shape),
      _direction_coord(other._direction_coord),
      _grid_frame(other._grid_frame),
      _map_frame(other._map_frame),
      _needs_rotation(other._needs_rotation),
      _rotation(other._rotation) {}

void CartaHealpixImage::ReadMap() {
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
    _maps.clear();
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
        _maps.push_back(map);
    }

    fits_close_file(fptr, &status);
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
        bool zenithal = IsZenithal(proj_type);
        width = _options.width > 0 ? _options.width : static_cast<int>(std::ceil((zenithal ? 90.0 : 360.0) / scale));
        height = _options.height > 0 ? _options.height : static_cast<int>(std::ceil((zenithal ? 90.0 : 180.0) / scale));

        double lon0 = std::isnan(_options.center_lon_deg) ? 0.0 : _options.center_lon_deg;
        double lat0 = std::isnan(_options.center_lat_deg) ? 0.0 : _options.center_lat_deg;

        casacore::Matrix<double> xform(2, 2, 0.0);
        xform.diagonal() = 1.0;
        // Reference pixel at the image centre (0-based), longitude increasing to the left
        _direction_coord = casacore::DirectionCoordinate(_grid_frame, casacore::Projection(proj_type), lon0 * DEG2RAD, lat0 * DEG2RAD,
            -scale * DEG2RAD, scale * DEG2RAD, xform, (width - 1) / 2.0, (height - 1) / 2.0);
        coord_sys.addCoordinate(_direction_coord);
    }

    // Column axis
    const size_t ncol = _maps.size();
    if (ncol > 1) {
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
    return !_maps.empty();
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
    const size_t npoints = nx * ny;

    if (!buffer.shape().isEqual(length)) {
        buffer.resize(length);
    }
    if (npoints == 0) {
        return false;
    }

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

    // The section is processed in chunks of points, in parallel, so that a full-resolution plane of a large grid
    // (tens of millions of pixels when matched to a survey image) never needs its coordinates held in memory at
    // once. Each chunk converts with its own copy of the direction coordinate (casacore's is not thread-safe).
    const size_t chunk_size(1 << 16);
    const size_t block_size(4096); // casacore's DirectionCoordinate::toWorldMany returns garbage for the points after
                                   // the first off-sky pixel of a batch, so any block with a failure is redone point by point
    const size_t num_chunks = (npoints + chunk_size - 1) / chunk_size;

    bool delete_it;
    float* data = buffer.getStorage(delete_it);

#pragma omp parallel for schedule(dynamic)
    for (size_t chunk = 0; chunk < num_chunks; ++chunk) {
        const size_t chunk_start = chunk * chunk_size;
        const size_t nchunk = std::min(chunk_size, npoints - chunk_start);
        casacore::DirectionCoordinate direction_copy(direction);
        casacore::Matrix<casacore::Double> pixel(2, nchunk, 0.0);
        std::vector<double> lon(nchunk), lat(nchunk);
        std::vector<bool> valid(nchunk, false);
        casacore::Vector<casacore::Double> one_pixel(2), one_world(2);

        // Grid pixel coordinates (0-based) of the chunk's points, x fastest
        for (size_t k = 0; k < nchunk; ++k) {
            const size_t p = chunk_start + k;
            pixel(lon_row, k) = blc(0) + (p % nx) * inc(0);
            pixel(lat_row, k) = blc(1) + (p / nx) * inc(1);
        }

        // Convert to world (radians)
        for (size_t block_start = 0; block_start < nchunk; block_start += block_size) {
            const size_t n = std::min(nchunk, block_start + block_size) - block_start;
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

        for (size_t k = 0; k < nchunk; ++k) {
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
            for (size_t k = 0; k < nchunk; ++k) {
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

        for (size_t c = 0; c < ncol; ++c) {
            size_t column = length.size() > 2 ? blc(2) + c * inc(2) : 0;
            Lookup(lon.data(), lat.data(), valid, nchunk, column, data + c * npoints + chunk_start);
        }
    }

    buffer.putStorage(data, delete_it);
    // As for CartaFitsImage, the return value is used by FileLoader::GetSlice as a success flag
    return true;
}

void CartaHealpixImage::Lookup(
    const double* lon, const double* lat, const std::vector<bool>& valid, size_t npoints, size_t column, float* out) const {
    if (column >= _maps.size()) {
        std::fill(out, out + npoints, FLOAT_NAN);
        return;
    }
    const std::vector<float>& map = *_maps[column];

    for (size_t p = 0; p < npoints; ++p) {
        if (!valid[p]) {
            out[p] = FLOAT_NAN;
            continue;
        }
        double theta = PI / 2.0 - lat[p];
        double phi = lon[p];
        if (phi < 0.0) {
            phi += 2.0 * PI;
        }
        pointing ptg(theta, phi);

        if (_options.bilinear) {
            fix_arr<int64_t, 4> pix;
            fix_arr<double, 4> wgt;
            _healpix->get_interpol(ptg, pix, wgt);
            double sum(0.0), wsum(0.0);
            for (int n = 0; n < 4; ++n) {
                float v = map[pix[n]];
                if (std::isfinite(v)) {
                    sum += wgt[n] * v;
                    wsum += wgt[n];
                }
            }
            out[p] = wsum > 0.0 ? static_cast<float>(sum / wsum) : FLOAT_NAN;
        } else {
            out[p] = map[_healpix->ang2pix(ptg)];
        }
    }
}
