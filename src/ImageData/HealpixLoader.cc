/* This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-backend
   Copyright 2018- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
   Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
   SPDX-License-Identifier: GPL-3.0-or-later
*/

#include "HealpixLoader.h"

#include <algorithm>
#include <cmath>
#include <limits>

#include "Logger/Logger.h"
#include "Timer/Timer.h"
#include "Util/Image.h"
#include "Util/Stokes.h"

namespace carta {

HealpixLoader::HealpixLoader(const std::string& filename) : FileLoader(filename), _hdu_num(0), _options_changed(false) {}

bool HealpixLoader::CanOpenFile(std::string& error) {
    if (!CartaHealpixImage::IsHealpixFile(_filename)) {
        error = "Not a HEALPix map.";
        return false;
    }
    return true;
}

void HealpixLoader::SetProjectionOptions(const HealpixProjectionOptions& options) {
    _options = options;
    _options_changed = true;
}

bool HealpixLoader::GetMapInfo(const std::string& hdu, HealpixMapInfo& info, std::string& error) {
    return CartaHealpixImage::ReadMapInfo(_filename, ResolveHdu(hdu), info, error);
}

unsigned int HealpixLoader::ResolveHdu(const std::string& hdu) {
    // HEALPix maps live in an extension; an empty or primary hdu means the first HEALPix extension
    unsigned int hdu_num = FileInfo::GetFitsHdu(hdu);
    if (hdu_num == 0) {
        unsigned int first_hdu(1);
        CartaHealpixImage::IsHealpixFile(_filename, &first_hdu);
        hdu_num = first_hdu;
    }
    return hdu_num;
}

void HealpixLoader::AllocateImage(const std::string& hdu) {
    unsigned int hdu_num = ResolveHdu(hdu);

    if (!_image || (hdu_num != _hdu_num) || _options_changed) {
        Timer t;
        _image.reset(new CartaHealpixImage(_filename, hdu_num, _options));
        _options_changed = false;

        if (!_image) {
            throw(casacore::AipsError("Error loading HEALPix map."));
        }

        _hdu = hdu;
        _hdu_num = hdu_num;
        _image_shape = _image->shape();
        _num_dims = _image_shape.size();
        _has_pixel_mask = false;
        _coord_sys = std::shared_ptr<casacore::CoordinateSystem>(static_cast<casacore::CoordinateSystem*>(_image->coordinates().clone()));
        _data_type = casacore::DataType::TpFloat;

        auto healpix_image = dynamic_cast<CartaHealpixImage*>(_image.get());
        const auto& info = healpix_image->MapInfo();
        spdlog::info("HEALPix map {}: NSIDE={} {} {} ({} {}{}), grid {}x{} in {:.1f} ms", _filename, info.nside,
            info.nested ? "NESTED" : "RING", info.coordsys.empty() ? "?" : info.coordsys, healpix_image->NumMaps(),
            info.image_format ? "map" : "column", healpix_image->NumMaps() == 1 ? "" : "s", _image_shape(0), _image_shape(1),
            t.Elapsed().ms());
    }
}

CartaHealpixImage* HealpixLoader::HealpixImage() const {
    return IsTransposed() ? nullptr : dynamic_cast<CartaHealpixImage*>(_image.get());
}

int HealpixLoader::MapIndex(int z, int stokes) const {
    // The third image axis (if any) is the map axis: Stokes for I/Q/U/V columns, else the spectral/linear/column axis
    if (_num_dims < 3) {
        return 0;
    }
    if (_axes.stokes == 2) {
        return stokes;
    }
    return z;
}

bool HealpixLoader::UseTileCache() const {
    return HealpixImage() != nullptr;
}

bool HealpixLoader::HasMip(int mip) const {
    return (mip >= 2) && (HealpixImage() != nullptr);
}

bool HealpixLoader::GetDownsampledRasterData(
    std::vector<float>& data, int z, int stokes, CARTA::ImageBounds& bounds, int mip, std::mutex& image_mutex) {
    auto image = HealpixImage();
    if (!image || mip < 1 || Stokes::IsComputed(stokes)) {
        return false;
    }
    // Downsampled pixel i covers full-resolution pixels [i mip, (i + 1) mip): sample the map at the block centre
    const int xmin = std::ceil((float)bounds.x_min() / mip);
    const int ymin = std::ceil((float)bounds.y_min() / mip);
    const int xmax = std::ceil((float)bounds.x_max() / mip);
    const int ymax = std::ceil((float)bounds.y_max() / mip);
    const int w = xmax - xmin;
    const int h = ymax - ymin;
    if (w <= 0 || h <= 0) {
        return false;
    }
    data.resize(static_cast<size_t>(w) * h);
    const double offset = (mip - 1) / 2.0;
    image->Sample(xmin * mip + offset, mip, w, ymin * mip + offset, mip, h, MapIndex(z, stokes), 1, 1, data.data());
    return true;
}

bool HealpixLoader::GetChunk(
    std::vector<float>& data, int& data_width, int& data_height, int min_x, int min_y, int z, int stokes, std::mutex& image_mutex) {
    if (!HealpixImage()) {
        return false;
    }
    data_width = std::min(CHUNK_SIZE, (int)_dims.width - min_x);
    data_height = std::min(CHUNK_SIZE, (int)_dims.height - min_y);
    if (data_width <= 0 || data_height <= 0) {
        return false;
    }

    // Through GetSlice, so that computed Stokes images (e.g. polarised intensity from Q and U) work as for other loaders
    StokesSource stokes_source(stokes, AxisRange(z), AxisRange(min_x, min_x + data_width - 1), AxisRange(min_y, min_y + data_height - 1));
    if (!stokes_source.IsOriginalImage()) {
        stokes = 0;
        z = 0;
        min_x = 0;
        min_y = 0;
    }
    casacore::Slicer slicer;
    if (_num_dims == 3) {
        slicer = casacore::Slicer(
            casacore::IPosition(3, min_x, min_y, _axes.stokes == 2 ? stokes : z), casacore::IPosition(3, data_width, data_height, 1));
    } else {
        slicer = casacore::Slicer(casacore::IPosition(2, min_x, min_y), casacore::IPosition(2, data_width, data_height));
    }
    data.resize(static_cast<size_t>(data_width) * data_height);
    casacore::Array<float> tmp(slicer.length(), data.data(), casacore::StorageInitPolicy::SHARE);
    try {
        GetSlice(tmp, StokesSlicer(stokes_source, slicer));
    } catch (casacore::AipsError& err) {
        spdlog::error("HEALPix chunk read failed: {}", err.getMesg());
        return false;
    }
    return true;
}

void HealpixLoader::LoadImageStats(bool) {
    // Size the per-plane entries; each is computed from the map on first use (GetImageStats)
    _z_stats.assign(_dims.num_stokes, std::vector<FileInfo::ImageStats>(_dims.depth));
    _cube_stats.assign(_dims.num_stokes, FileInfo::ImageStats());
}

FileInfo::ImageStats& HealpixLoader::GetImageStats(int current_stokes, int z) {
    auto image = HealpixImage();
    if (!image || Stokes::IsComputed(current_stokes) || z < 0 || current_stokes < 0 || current_stokes >= (int)_z_stats.size() ||
        z >= (int)_z_stats[current_stokes].size()) {
        return FileLoader::GetImageStats(current_stokes, z);
    }
    std::lock_guard<std::mutex> guard(_stats_mutex);
    auto& stats = _z_stats[current_stokes][z];
    if (!stats.valid) {
        Timer t;
        ComputeMapStats(MapIndex(z, current_stokes), stats);
        spdlog::performance("HEALPix map statistics for plane {} in {:.3f} ms", MapIndex(z, current_stokes), t.Elapsed().ms());
    }
    return stats;
}

void HealpixLoader::ComputeMapStats(size_t map, FileInfo::ImageStats& stats) const {
    // Statistics of the map's own pixels (every HEALPix pixel once, so independent of the projection and grid size)
    auto image = HealpixImage();
    const size_t npix = image->Npix();
    auto held = image->MapValues(map);
    auto value_at = [&](size_t pix) { return held ? (*held)[pix] : image->MapValue(map, pix); };

    double sum(0.0), sum_sq(0.0);
    float min_val(std::numeric_limits<float>::max()), max_val(std::numeric_limits<float>::lowest());
    size_t num_pixels(0);
    for (size_t pix = 0; pix < npix; ++pix) {
        float v = value_at(pix);
        if (!std::isfinite(v)) {
            continue;
        }
        ++num_pixels;
        sum += v;
        sum_sq += static_cast<double>(v) * v;
        min_val = std::min(min_val, v);
        max_val = std::max(max_val, v);
    }

    auto& basic = stats.basic_stats;
    basic[CARTA::StatsType::NumPixels] = num_pixels;
    basic[CARTA::StatsType::NanCount] = npix - num_pixels;
    if (num_pixels == 0) {
        basic[CARTA::StatsType::Min] = basic[CARTA::StatsType::Max] = std::nan("");
        stats.histogram_bins.clear();
        stats.valid = true;
        stats.full = true;
        return;
    }
    basic[CARTA::StatsType::Sum] = sum;
    basic[CARTA::StatsType::SumSq] = sum_sq;
    basic[CARTA::StatsType::Min] = min_val;
    basic[CARTA::StatsType::Max] = max_val;
    basic[CARTA::StatsType::Mean] = sum / num_pixels;
    basic[CARTA::StatsType::Sigma] = num_pixels > 1 ? std::sqrt(std::max(0.0, (sum_sq - sum * sum / num_pixels) / (num_pixels - 1))) : 0.0;
    basic[CARTA::StatsType::RMS] = std::sqrt(sum_sq / num_pixels);
    basic[CARTA::StatsType::Extrema] = std::abs(min_val) > std::abs(max_val) ? min_val : max_val;

    // Histogram over [min, max], as many bins as the automatic choice for an image of this many pixels (capped)
    const int num_bins = std::max(2, std::min(2048, static_cast<int>(std::sqrt(static_cast<double>(num_pixels)))));
    stats.histogram_bins.assign(num_bins, 0);
    const double width = (static_cast<double>(max_val) - min_val) / num_bins;
    for (size_t pix = 0; pix < npix; ++pix) {
        float v = value_at(pix);
        if (!std::isfinite(v)) {
            continue;
        }
        int bin = width > 0 ? static_cast<int>((v - min_val) / width) : 0;
        stats.histogram_bins[std::min(num_bins - 1, std::max(0, bin))]++;
    }
    stats.valid = true;
    stats.full = true;
}

} // namespace carta
