/* This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-backend
   Copyright 2018- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
   Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
   SPDX-License-Identifier: GPL-3.0-or-later
*/

#include "HealpixLoader.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>

#include "Logger/Logger.h"
#include "Timer/Timer.h"
#include "Util/Image.h"
#include "Util/Nan.h"
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
    // Statistics of the map's own pixels (every HEALPix pixel once, so independent of the projection and grid size),
    // read in blocks
    auto image = HealpixImage();
    const size_t npix = image->Npix();
    const size_t block_size(1 << 18);
    std::vector<float> block(std::min(block_size, npix));
    auto for_each_value = [&](const std::function<void(float)>& f) {
        for (size_t start = 0; start < npix; start += block_size) {
            const size_t n = std::min(block_size, npix - start);
            image->ReadMapRange(map, start, n, block.data());
            for (size_t i = 0; i < n; ++i) {
                if (std::isfinite(block[i])) {
                    f(block[i]);
                }
            }
        }
    };

    double sum(0.0), sum_sq(0.0);
    float min_val(std::numeric_limits<float>::max()), max_val(std::numeric_limits<float>::lowest());
    size_t num_pixels(0);
    for_each_value([&](float v) {
        ++num_pixels;
        sum += v;
        sum_sq += static_cast<double>(v) * v;
        min_val = std::min(min_val, v);
        max_val = std::max(max_val, v);
    });

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
    for_each_value([&](float v) {
        int bin = width > 0 ? static_cast<int>((v - min_val) / width) : 0;
        stats.histogram_bins[std::min(num_bins - 1, std::max(0, bin))]++;
    });
    stats.valid = true;
    stats.full = true;
}

bool HealpixLoader::HasSpectralMaps() const {
    // A cube whose third axis (the maps) is the spectral axis of the frame, not a Stokes axis
    return HealpixImage() != nullptr && _num_dims == 3 && _axes.z == 2 && _dims.depth > 1;
}

bool HealpixLoader::GetCursorSpectralData(
    std::vector<float>& data, int stokes, int cursor_x, int count_x, int cursor_y, int count_y, std::mutex& image_mutex) {
    auto image = HealpixImage();
    if (!HasSpectralMaps() || Stokes::IsComputed(stokes) || count_x < 1 || count_y < 1) {
        return false;
    }
    if (cursor_x < 0 || cursor_y < 0 || cursor_x + count_x > (int)_dims.width || cursor_y + count_y > (int)_dims.height) {
        return false;
    }

    // All maps at the pixels, then reordered to z fastest (then y, then x) as the swizzled HDF5 data is laid out
    const size_t depth = _dims.depth;
    const size_t npoints = static_cast<size_t>(count_x) * count_y;
    std::vector<float> sampled(npoints * depth);
    std::lock_guard<std::mutex> guard(image_mutex);
    try {
        image->Sample(cursor_x, 1, count_x, cursor_y, 1, count_y, 0, 1, depth, sampled.data());
    } catch (casacore::AipsError& err) {
        spdlog::error("HEALPix cursor spectral profile failed: {}", err.getMesg());
        return false;
    }
    data.resize(npoints * depth);
    for (int x = 0; x < count_x; ++x) {
        for (int y = 0; y < count_y; ++y) {
            const size_t point = y * count_x + x;
            float* profile = data.data() + (static_cast<size_t>(x) * count_y + y) * depth;
            for (size_t z = 0; z < depth; ++z) {
                profile[z] = sampled[z * npoints + point];
            }
        }
    }
    return true;
}

bool HealpixLoader::UseRegionSpectralData(const casacore::IPosition&, std::mutex&) {
    return HasSpectralMaps();
}

bool HealpixLoader::GetRegionSpectralData(int region_id, const AxisRange& spectral_range, int stokes,
    const casacore::ArrayLattice<casacore::Bool>& mask, const casacore::IPosition& origin, std::mutex& image_mutex,
    std::map<CARTA::StatsType, std::vector<double>>& results, float& progress) {
    // Per-channel statistics of the region, accumulated a band of rows at a time: each call converts the band's grid
    // points to map pixels once, reads them from every channel's map, and returns the stats so far with the progress.
    // The calling function keeps calling until the progress is complete.
    auto image = HealpixImage();
    if (!HasSpectralMaps() || Stokes::IsComputed(stokes) || mask.shape().size() < 2 || origin.size() < 2) {
        return false;
    }

    AxisRange z_range(spectral_range.from, spectral_range.to);
    if (z_range.to == ALL_Z || z_range.to >= (int)_dims.depth) {
        z_range.to = _dims.depth - 1;
    }
    if (z_range.from < 0 || z_range.from > z_range.to) {
        return false;
    }
    const int depth = z_range.to - z_range.from + 1;
    const casacore::IPosition mask_shape(mask.shape());
    const int width = mask_shape(0);
    const int height = mask_shape(1);
    const casacore::Array<casacore::Bool> mask_array = mask.get();
    bool delete_mask;
    const casacore::Bool* mask_data = mask_array.getStorage(delete_mask);
    const size_t mask_size = static_cast<size_t>(width) * height;

    // A cheap signature of the mask: which pixels of the bounding box belong to the region
    size_t signature(mask_size);
    for (size_t i = 0; i < mask_size; ++i) {
        if (mask_data[i]) {
            signature = signature * 1099511628211ULL + i;
        }
    }

    std::lock_guard<std::mutex> guard(_spectral_mutex);
    auto id = FileInfo::RegionStatsId(region_id, stokes);
    auto& state = _region_spectral[id];
    if (!state.origin.isEqual(origin) || !state.shape.isEqual(mask_shape) || state.mask_signature != signature ||
        state.z_from != z_range.from || state.z_to != z_range.to) {
        state = RegionSpectralProgress();
        state.origin = origin;
        state.shape = mask_shape;
        state.mask_signature = signature;
        state.z_from = z_range.from;
        state.z_to = z_range.to;
        state.num_pixels.assign(depth, 0);
        state.sum.assign(depth, 0.0);
        state.sum_sq.assign(depth, 0.0);
        state.min.assign(depth, std::numeric_limits<double>::max());
        state.max.assign(depth, std::numeric_limits<double>::lowest());
    }

    if (!state.completed) {
        // The next band of rows: enough grid points to be worth a read from every map, few enough to answer promptly
        const size_t band_points(1 << 15);
        const size_t rows = std::max<size_t>(1, band_points / std::max(1, width));
        const size_t row_end = std::min<size_t>(height, state.next_row + rows);

        // The band's region points (0-based grid pixel coordinates), converted to map pixels once
        std::vector<double> x, y;
        for (size_t row = state.next_row; row < row_end; ++row) {
            for (int col = 0; col < width; ++col) {
                if (mask_data[row * width + col]) {
                    x.push_back(origin(0) + col);
                    y.push_back(origin(1) + row);
                }
            }
        }
        if (!x.empty()) {
            HealpixLookupPlan plan;
            std::vector<float> values(x.size());
            std::lock_guard<std::mutex> image_guard(image_mutex);
            try {
                image->PlanLookup(x.data(), y.data(), x.size(), plan);
                for (int z = z_range.from; z <= z_range.to; ++z) {
                    image->LookupPlanned(plan, MapIndex(z, stokes), values.data());
                    const int zi = z - z_range.from;
                    for (float v : values) {
                        if (std::isfinite(v)) {
                            state.num_pixels[zi]++;
                            state.sum[zi] += v;
                            state.sum_sq[zi] += static_cast<double>(v) * v;
                            state.min[zi] = std::min<double>(state.min[zi], v);
                            state.max[zi] = std::max<double>(state.max[zi], v);
                        }
                    }
                }
            } catch (casacore::AipsError& err) {
                spdlog::error("HEALPix region spectral profile failed: {}", err.getMesg());
                mask_array.freeStorage(mask_data, delete_mask);
                return false;
            }
        }
        state.next_row = row_end;
        state.completed = (row_end >= static_cast<size_t>(height));

        // Stats from the accumulators so far
        const double beam_area = CalculateBeamArea();
        const bool has_flux = std::isfinite(beam_area) && beam_area > 0;
        auto& stats = state.stats;
        for (auto type : {CARTA::StatsType::NumPixels, CARTA::StatsType::NanCount, CARTA::StatsType::Sum, CARTA::StatsType::Mean,
                 CARTA::StatsType::RMS, CARTA::StatsType::Sigma, CARTA::StatsType::SumSq, CARTA::StatsType::Min, CARTA::StatsType::Max,
                 CARTA::StatsType::Extrema}) {
            stats[type].assign(depth, DOUBLE_NAN);
        }
        if (has_flux) {
            stats[CARTA::StatsType::FluxDensity].assign(depth, DOUBLE_NAN);
        }
        for (int zi = 0; zi < depth; ++zi) {
            const double n = state.num_pixels[zi];
            stats[CARTA::StatsType::NumPixels][zi] = n;
            stats[CARTA::StatsType::NanCount][zi] = 0; // off-sky and unseen points are not counted either way
            if (n == 0) {
                continue;
            }
            const double sum = state.sum[zi], sum_sq = state.sum_sq[zi];
            stats[CARTA::StatsType::Sum][zi] = sum;
            stats[CARTA::StatsType::SumSq][zi] = sum_sq;
            stats[CARTA::StatsType::Mean][zi] = sum / n;
            stats[CARTA::StatsType::RMS][zi] = std::sqrt(sum_sq / n);
            stats[CARTA::StatsType::Sigma][zi] = n > 1 ? std::sqrt(std::max(0.0, (sum_sq - sum * sum / n) / (n - 1))) : 0.0;
            stats[CARTA::StatsType::Min][zi] = state.min[zi];
            stats[CARTA::StatsType::Max][zi] = state.max[zi];
            stats[CARTA::StatsType::Extrema][zi] = std::abs(state.min[zi]) > std::abs(state.max[zi]) ? state.min[zi] : state.max[zi];
            if (has_flux) {
                stats[CARTA::StatsType::FluxDensity][zi] = sum / beam_area;
            }
        }
    }
    mask_array.freeStorage(mask_data, delete_mask);

    results = state.stats;
    progress = state.completed ? 1.0f : static_cast<float>(state.next_row) / height;
    if (state.completed && region_id <= TEMP_REGION_ID) {
        _region_spectral.erase(id); // temporary regions (spatial profiles, PV images) are not kept
    }
    return true;
}

} // namespace carta
