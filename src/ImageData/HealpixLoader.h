/* This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-backend
   Copyright 2018- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
   Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
   SPDX-License-Identifier: GPL-3.0-or-later
*/

//# HealpixLoader.h : loader for HEALPix maps (FITS binary tables), reprojected on the fly by CartaHealpixImage

#ifndef CARTA_SRC_IMAGEDATA_HEALPIXLOADER_H_
#define CARTA_SRC_IMAGEDATA_HEALPIXLOADER_H_

#include <mutex>

#include "CartaHealpixImage.h"
#include "FileLoader.h"

namespace carta {

class HealpixLoader : public FileLoader {
public:
    HealpixLoader(const std::string& filename);

    bool CanOpenFile(std::string& error) override;

    /**
     * @brief Set the projection/grid the map is reprojected onto. Takes effect on the next OpenFile;
     * an already-open image with different options is discarded so it is rebuilt.
     */
    void SetProjectionOptions(const HealpixProjectionOptions& options);
    const HealpixProjectionOptions& GetProjectionOptions() const {
        return _options;
    }

    /**
     * @brief Metadata of the HEALPix HDU (without reading the map), for file info.
     */
    bool GetMapInfo(const std::string& hdu, HealpixMapInfo& info, std::string& error);

    // Tiles and downsampled tiles are reprojected on demand at any level, so the frame never needs a full-resolution
    // plane of the grid (which may be extended far beyond the matched image). Not for a transposed view.
    bool UseTileCache() const override;
    bool HasMip(int mip) const override;
    bool GetDownsampledRasterData(
        std::vector<float>& data, int z, int stokes, CARTA::ImageBounds& bounds, int mip, std::mutex& image_mutex) override;
    bool GetChunk(std::vector<float>& data, int& data_width, int& data_height, int min_x, int min_y, int z, int stokes,
        std::mutex& image_mutex) override;

    // Per-plane statistics and histogram come from the HEALPix map's own pixels (the whole sky, not the grid)
    void LoadImageStats(bool load_percentiles = false) override;
    FileInfo::ImageStats& GetImageStats(int current_stokes, int z) override;

private:
    void AllocateImage(const std::string& hdu) override;
    unsigned int ResolveHdu(const std::string& hdu);
    CartaHealpixImage* HealpixImage() const;
    int MapIndex(int z, int stokes) const;
    void ComputeMapStats(size_t map, FileInfo::ImageStats& stats) const;

    HealpixProjectionOptions _options;
    unsigned int _hdu_num;
    bool _options_changed;
    std::mutex _stats_mutex;
};

} // namespace carta

#endif // CARTA_SRC_IMAGEDATA_HEALPIXLOADER_H_
