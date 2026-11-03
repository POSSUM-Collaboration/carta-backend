/* This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-backend
   Copyright 2018- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
   Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
   SPDX-License-Identifier: GPL-3.0-or-later
*/

//# TransposedImage.h : a lazily axis-permuted view of another casacore::ImageInterface, so that any file axis can
//# be rendered or animated without rewriting the image

#ifndef CARTA_SRC_IMAGEDATA_TRANSPOSEDIMAGE_H_
#define CARTA_SRC_IMAGEDATA_TRANSPOSEDIMAGE_H_

#include <memory>
#include <string>
#include <vector>

#include <casacore/casa/Arrays/IPosition.h>
#include <casacore/casa/Arrays/Vector.h>
#include <casacore/casa/BasicSL/String.h>
#include <casacore/images/Images/ImageInterface.h>

namespace carta {

/**
 * @brief Presents another image with its axes permuted: image axis i is file axis order[i].
 *
 * Slices are read from the wrapped image with the slicer mapped back to file order and the data reordered,
 * the coordinate system is transposed accordingly, and units, image info and misc info are copied. The view
 * is read-only and holds a shared pointer to the wrapped image.
 */
class TransposedImage : public casacore::ImageInterface<float> {
public:
    TransposedImage(std::shared_ptr<casacore::ImageInterface<float>> image, const std::vector<int>& axis_order);
    TransposedImage(const TransposedImage& other);
    ~TransposedImage() override = default;

    /// Whether order is a permutation of 0..ndim-1
    static bool IsValidOrder(const std::vector<int>& axis_order, size_t ndim);
    /// Whether order leaves every axis in place (or is empty)
    static bool IsIdentityOrder(const std::vector<int>& axis_order);

    /**
     * @brief Renumber a FITS axis keyword (CTYPEn, CRVALn, NAXISn, PCi_j, CDi_j, PVi_m, ... with optional alternate
     * WCS letter) from file axis numbers to image axis numbers. inverse_order[f] is the 0-based image axis of file axis f.
     * Keywords that do not refer to an axis are returned unchanged.
     */
    static std::string RenumberAxisKeyword(const std::string& keyword, const std::vector<int>& inverse_order);
    /// Renumber the keywords of 80-character FITS header cards in place
    static void RenumberHeaderCards(casacore::Vector<casacore::String>& headers, const std::vector<int>& inverse_order);

    std::shared_ptr<casacore::ImageInterface<float>> Inner() const {
        return _image;
    }
    const std::vector<int>& AxisOrder() const {
        return _order; // image axis i = file axis _order[i]
    }
    const std::vector<int>& InverseOrder() const {
        return _inverse; // file axis f = image axis _inverse[f]
    }

    // casacore::ImageInterface
    casacore::String imageType() const override;
    casacore::String name(bool stripPath = false) const override;
    casacore::IPosition shape() const override;
    casacore::Bool ok() const override;
    casacore::DataType dataType() const override;
    casacore::Bool doGetSlice(casacore::Array<float>& buffer, const casacore::Slicer& section) override;
    void doPutSlice(const casacore::Array<float>& buffer, const casacore::IPosition& where, const casacore::IPosition& stride) override;
    const casacore::LatticeRegion* getRegionPtr() const override;
    casacore::ImageInterface<float>* cloneII() const override;
    void resize(const casacore::TiledShape& newShape) override;
    casacore::uInt advisedMaxPixels() const override;
    casacore::IPosition doNiceCursorShape(casacore::uInt maxPixels) const override;
    casacore::Bool isMasked() const override;
    casacore::Bool hasPixelMask() const override;
    casacore::Bool doGetMaskSlice(casacore::Array<bool>& buffer, const casacore::Slicer& section) override;
    void tempClose() override;
    void reopen() override;

private:
    casacore::Slicer ToFileOrder(const casacore::Slicer& section) const;
    casacore::IPosition ToImageOrder(const casacore::IPosition& file_position) const;

    std::shared_ptr<casacore::ImageInterface<float>> _image;
    std::vector<int> _order;
    std::vector<int> _inverse;
    casacore::IPosition _shape;
};

} // namespace carta

#endif // CARTA_SRC_IMAGEDATA_TRANSPOSEDIMAGE_H_
