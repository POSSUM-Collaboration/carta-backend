/* This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-backend
   Copyright 2018- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
   Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
   SPDX-License-Identifier: GPL-3.0-or-later
*/

#include "TransposedImage.h"

#include <algorithm>
#include <cctype>
#include <numeric>

#include <casacore/casa/Arrays/ArrayUtil.h>
#include <casacore/casa/Exceptions/Error.h>
#include <casacore/coordinates/Coordinates/CoordinateSystem.h>

using namespace carta;

bool TransposedImage::IsValidOrder(const std::vector<int>& axis_order, size_t ndim) {
    if (axis_order.size() != ndim) {
        return false;
    }
    std::vector<bool> seen(ndim, false);
    for (auto axis : axis_order) {
        if ((axis < 0) || (axis >= static_cast<int>(ndim)) || seen[axis]) {
            return false;
        }
        seen[axis] = true;
    }
    return true;
}

bool TransposedImage::IsIdentityOrder(const std::vector<int>& axis_order) {
    for (size_t i = 0; i < axis_order.size(); ++i) {
        if (axis_order[i] != static_cast<int>(i)) {
            return false;
        }
    }
    return true;
}

std::string TransposedImage::RenumberAxisKeyword(const std::string& keyword, const std::vector<int>& inverse_order) {
    // Axis keywords: PREFIX + axis number [+ alternate letter], or PREFIX + i + '_' + j (both axes for PC/CD, first for PV/PS)
    static const std::vector<std::string> single_axis = {
        "CTYPE", "CRVAL", "CDELT", "CRPIX", "CUNIT", "CROTA", "CNAME", "CRDER", "CSYER", "NAXIS"};
    static const std::vector<std::string> matrix_both = {"PC", "CD"};
    static const std::vector<std::string> matrix_first = {"PV", "PS"};

    auto map_axis = [&](const std::string& digits, std::string& out) {
        int file_axis = std::stoi(digits) - 1;
        if ((file_axis < 0) || (file_axis >= static_cast<int>(inverse_order.size()))) {
            return false;
        }
        out = std::to_string(inverse_order[file_axis] + 1);
        return true;
    };
    auto split_suffix = [](const std::string& rest, std::string& digits, std::string& alt) {
        size_t i = 0;
        while (i < rest.size() && std::isdigit(static_cast<unsigned char>(rest[i]))) {
            ++i;
        }
        digits = rest.substr(0, i);
        alt = rest.substr(i);
        return !digits.empty() && (alt.empty() || (alt.size() == 1 && std::isupper(static_cast<unsigned char>(alt[0]))));
    };

    for (const auto& prefix : single_axis) {
        if (keyword.rfind(prefix, 0) == 0 && keyword.size() > prefix.size()) {
            std::string digits, alt, mapped;
            if (split_suffix(keyword.substr(prefix.size()), digits, alt) && map_axis(digits, mapped)) {
                return prefix + mapped + alt;
            }
            return keyword;
        }
    }
    for (int which = 0; which < 2; ++which) {
        const auto& prefixes = which == 0 ? matrix_both : matrix_first;
        for (const auto& prefix : prefixes) {
            if (keyword.rfind(prefix, 0) == 0 && keyword.size() > prefix.size()) {
                auto rest = keyword.substr(prefix.size());
                auto underscore = rest.find('_');
                if (underscore == std::string::npos) {
                    return keyword;
                }
                std::string first = rest.substr(0, underscore), second = rest.substr(underscore + 1);
                std::string digits2, alt, mapped1, mapped2;
                if (!split_suffix(second, digits2, alt) || !map_axis(first, mapped1)) {
                    return keyword;
                }
                if (which == 0) {
                    if (!map_axis(digits2, mapped2)) {
                        return keyword;
                    }
                } else {
                    mapped2 = digits2;
                }
                return prefix + mapped1 + "_" + mapped2 + alt;
            }
        }
    }
    return keyword;
}

void TransposedImage::RenumberHeaderCards(casacore::Vector<casacore::String>& headers, const std::vector<int>& inverse_order) {
    for (size_t i = 0; i < headers.size(); ++i) {
        std::string card(headers(i));
        if (card.size() < 9) {
            continue;
        }
        std::string keyword = card.substr(0, 8);
        keyword.erase(keyword.find_last_not_of(' ') + 1);
        std::string renumbered = RenumberAxisKeyword(keyword, inverse_order);
        if (renumbered != keyword) {
            renumbered.resize(8, ' ');
            headers(i) = casacore::String(renumbered + card.substr(8));
        }
    }
}

TransposedImage::TransposedImage(std::shared_ptr<casacore::ImageInterface<float>> image, const std::vector<int>& axis_order)
    : casacore::ImageInterface<float>(), _image(image), _order(axis_order) {
    if (!_image) {
        throw(casacore::AipsError("TransposedImage: no image."));
    }
    const size_t ndim = _image->ndim();
    if (!IsValidOrder(_order, ndim)) {
        throw(casacore::AipsError("Axis order must list each of the " + std::to_string(ndim) + " image axes exactly once (0-based)."));
    }
    _inverse.assign(ndim, 0);
    for (size_t i = 0; i < ndim; ++i) {
        _inverse[_order[i]] = i;
    }

    const auto file_shape = _image->shape();
    _shape.resize(ndim);
    for (size_t i = 0; i < ndim; ++i) {
        _shape(i) = file_shape(_order[i]);
    }

    // Transpose the coordinate system: new pixel axis i is old pixel axis order[i]; world axes follow when they
    // correspond one-to-one with pixel axes (always true for images read from files)
    casacore::CoordinateSystem coord_sys(_image->coordinates());
    casacore::Vector<casacore::Int> pixel_order(ndim), world_order(coord_sys.nWorldAxes());
    for (size_t i = 0; i < ndim; ++i) {
        pixel_order(i) = _order[i];
    }
    if (coord_sys.nWorldAxes() == ndim) {
        world_order = pixel_order;
    } else {
        std::iota(world_order.begin(), world_order.end(), 0);
    }
    coord_sys.transpose(world_order, pixel_order);
    setCoordinateInfo(coord_sys);

    setUnits(_image->units());
    setImageInfo(_image->imageInfo());
    setMiscInfo(_image->miscInfo());
}

TransposedImage::TransposedImage(const TransposedImage& other)
    : casacore::ImageInterface<float>(other), _image(other._image), _order(other._order), _inverse(other._inverse), _shape(other._shape) {}

casacore::Slicer TransposedImage::ToFileOrder(const casacore::Slicer& section) const {
    const size_t ndim = _order.size();
    casacore::IPosition start(ndim), end(ndim), stride(ndim);
    for (size_t i = 0; i < ndim; ++i) {
        start(_order[i]) = section.start()(i);
        end(_order[i]) = section.end()(i);
        stride(_order[i]) = section.stride()(i);
    }
    return casacore::Slicer(start, end, stride, casacore::Slicer::endIsLast);
}

casacore::IPosition TransposedImage::ToImageOrder(const casacore::IPosition& file_position) const {
    casacore::IPosition image_position(_order.size());
    for (size_t i = 0; i < _order.size(); ++i) {
        image_position(i) = file_position(_order[i]);
    }
    return image_position;
}

casacore::String TransposedImage::imageType() const {
    return "TransposedImage";
}

casacore::String TransposedImage::name(bool stripPath) const {
    return _image->name(stripPath);
}

casacore::IPosition TransposedImage::shape() const {
    return _shape;
}

casacore::Bool TransposedImage::ok() const {
    return _image->ok();
}

casacore::DataType TransposedImage::dataType() const {
    return casacore::TpFloat;
}

casacore::Bool TransposedImage::doGetSlice(casacore::Array<float>& buffer, const casacore::Slicer& section) {
    // Read the section in file order, then move the axes into image order
    casacore::Array<float> file_data;
    _image->getSlice(file_data, ToFileOrder(section));
    casacore::IPosition new_order(_order.size());
    for (size_t i = 0; i < _order.size(); ++i) {
        new_order(i) = _order[i];
    }
    casacore::Array<float> image_data = casacore::reorderArray(file_data, new_order, false);
    if (!buffer.shape().isEqual(image_data.shape())) {
        buffer.resize(image_data.shape());
    }
    buffer = image_data;
    return false;
}

casacore::Bool TransposedImage::doGetMaskSlice(casacore::Array<bool>& buffer, const casacore::Slicer& section) {
    if (!_image->isMasked()) {
        buffer.resize(section.length());
        buffer = true;
        return false;
    }
    casacore::Array<bool> file_mask;
    _image->getMaskSlice(file_mask, ToFileOrder(section));
    casacore::IPosition new_order(_order.size());
    for (size_t i = 0; i < _order.size(); ++i) {
        new_order(i) = _order[i];
    }
    casacore::Array<bool> image_mask = casacore::reorderArray(file_mask, new_order, false);
    if (!buffer.shape().isEqual(image_mask.shape())) {
        buffer.resize(image_mask.shape());
    }
    buffer = image_mask;
    return false;
}

void TransposedImage::doPutSlice(const casacore::Array<float>&, const casacore::IPosition&, const casacore::IPosition&) {
    throw(casacore::AipsError("TransposedImage is read-only."));
}

const casacore::LatticeRegion* TransposedImage::getRegionPtr() const {
    return nullptr;
}

casacore::ImageInterface<float>* TransposedImage::cloneII() const {
    return new TransposedImage(*this);
}

void TransposedImage::resize(const casacore::TiledShape&) {
    throw(casacore::AipsError("TransposedImage cannot be resized."));
}

casacore::uInt TransposedImage::advisedMaxPixels() const {
    return _image->advisedMaxPixels();
}

casacore::IPosition TransposedImage::doNiceCursorShape(casacore::uInt maxPixels) const {
    return ToImageOrder(_image->niceCursorShape(maxPixels));
}

casacore::Bool TransposedImage::isMasked() const {
    return _image->isMasked();
}

casacore::Bool TransposedImage::hasPixelMask() const {
    return false; // masks are served through doGetMaskSlice
}

void TransposedImage::tempClose() {
    _image->tempClose();
}

void TransposedImage::reopen() {
    _image->reopen();
}
