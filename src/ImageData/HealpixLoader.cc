/* This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-backend
   Copyright 2018- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
   Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
   SPDX-License-Identifier: GPL-3.0-or-later
*/

#include "HealpixLoader.h"

#include "Logger/Logger.h"
#include "Timer/Timer.h"

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

} // namespace carta
