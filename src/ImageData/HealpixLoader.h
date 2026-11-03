/* This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-backend
   Copyright 2018- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
   Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
   SPDX-License-Identifier: GPL-3.0-or-later
*/

//# HealpixLoader.h : loader for HEALPix maps (FITS binary tables), reprojected on the fly by CartaHealpixImage

#ifndef CARTA_SRC_IMAGEDATA_HEALPIXLOADER_H_
#define CARTA_SRC_IMAGEDATA_HEALPIXLOADER_H_

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

private:
    void AllocateImage(const std::string& hdu) override;
    unsigned int ResolveHdu(const std::string& hdu);

    HealpixProjectionOptions _options;
    unsigned int _hdu_num;
    bool _options_changed;
};

} // namespace carta

#endif // CARTA_SRC_IMAGEDATA_HEALPIXLOADER_H_
