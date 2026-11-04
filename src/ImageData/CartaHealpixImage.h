/* This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-backend
   Copyright 2018- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
   Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
   SPDX-License-Identifier: GPL-3.0-or-later
*/

//# CartaHealpixImage.h : HEALPix map (FITS binary table, or an image HDU holding one map per row) presented as a
//# casacore::ImageInterface by reprojecting it on the fly onto a regular pixel grid with a WCS of the caller's choosing.

#ifndef CARTA_SRC_IMAGEDATA_CARTAHEALPIXIMAGE_H_
#define CARTA_SRC_IMAGEDATA_CARTAHEALPIXIMAGE_H_

#include <memory>
#include <string>
#include <vector>

#include <casacore/casa/Arrays/Matrix.h>
#include <casacore/coordinates/Coordinates/CoordinateSystem.h>
#include <casacore/coordinates/Coordinates/DirectionCoordinate.h>
#include <casacore/images/Images/ImageInterface.h>
#include <casacore/measures/Measures/MDirection.h>

// healpix_cxx headers are only included in the .cc: its global typedefs clash with casacore's mirlib
template <typename I>
class T_Healpix_Base;

namespace carta {

class HealpixMapSource;

/**
 * @brief Metadata of a HEALPix HDU: a FITS binary table (one column per map), or an image HDU whose first axis is
 * the HEALPix pixel index and whose further axes enumerate maps (e.g. a spectral cube, one full-sky map per channel).
 */
struct HealpixMapInfo {
    long nside = 0;
    bool nested = false;                   ///< ORDERING = NESTED (else RING)
    bool explicit_index = false;           ///< INDXSCHM = EXPLICIT (partial sky, PIXEL column)
    std::string coordsys;                  ///< COORDSYS header: G, C/Q, E
    std::string object;                    ///< OBJECT header: FULLSKY or PARTIAL
    std::vector<std::string> columns;      ///< data column names (TTYPEn), excluding the PIXEL column; "MAP" for an image HDU
    std::vector<std::string> column_units; ///< TUNITn for the data columns (BUNIT for an image HDU)
    double bad_data = -1.6375e30;          ///< BAD_DATA header (HEALPix UNSEEN)

    // Image-format maps (PIXTYPE = HEALPIX on an image HDU with NAXIS1 = 12 NSIDE^2)
    bool image_format = false;
    long num_maps = 1;          ///< product of the axes after the first
    std::string map_axis_ctype; ///< WCS of the map axis (CTYPEn/CRVALn/CDELTn/CRPIXn/CUNITn, n >= 2), if any
    std::string map_axis_unit;
    double map_axis_crval = 0.0;
    double map_axis_cdelt = 1.0;
    double map_axis_crpix = 1.0;        ///< FITS 1-based
    std::string specsys;                ///< SPECSYS header, if any
    double rest_freq_hz = std::nan(""); ///< RESTFRQ / RESTFREQ, if any
    double bmaj_deg = std::nan("");     ///< BMAJ / BMIN / BPA, if any
    double bmin_deg = std::nan("");
    double bpa_deg = std::nan("");

    long npix() const;
    double resolution_deg() const; ///< mean pixel spacing: sqrt(4 pi / npix)
    casacore::MDirection::Types frame() const;
    bool has_map_axis() const {
        return image_format && num_maps > 1 && !map_axis_ctype.empty();
    }
};

/**
 * @brief Options describing the pixel grid and projection a HEALPix map is reprojected onto.
 *
 * Any field left at its default is filled in from the map: all-sky Mollweide in the map's own
 * coordinate frame, at the HEALPix pixel resolution. A matched coordinate system (from another
 * open image) overrides the projection fields entirely.
 */
struct HealpixProjectionOptions {
    std::string projection = "MOL";                              ///< WCS projection code (MOL, AIT, CAR, SFL, ZEA, TAN, SIN, HPX, ...)
    std::string frame;                                           ///< target celestial frame: G, C (ICRS/J2000), E; empty = map frame
    double center_lon_deg = std::nan("");                        ///< reference longitude; NaN = 0
    double center_lat_deg = std::nan("");                        ///< reference latitude; NaN = 0 (all-sky) or map centre
    double pixel_scale_deg = std::nan("");                       ///< NaN = HEALPix resolution
    int width = 0;                                               ///< 0 = enough to cover the sky (or the matched image's size)
    int height = 0;                                              ///< when matching, a size other than the image's keeps its grid centred
    bool bilinear = false;                                       ///< bilinear interpolation across the 4 nearest HEALPix pixels
    std::vector<std::string> columns;                            ///< empty = all data columns
    std::shared_ptr<casacore::CoordinateSystem> match_coord_sys; ///< reproject onto this system's direction coordinate
    casacore::IPosition match_shape;                             ///< and this pixel grid (render axes)
};

/**
 * @brief Read-only casacore image backed by a HEALPix map (a table held in memory, or an image HDU read lazily).
 *
 * Shape is [width, height] for a single map, or [width, height, nmaps] with a Stokes axis (when the table columns
 * are recognisable I/Q/U/V), a spectral or linear axis (image-format maps with a WCS on their map axis), or a
 * linear "COLUMN" axis otherwise. Pixel values
 * are produced in doGetSlice by converting each requested grid pixel to a sky direction with the
 * direction coordinate, rotating into the map's frame if necessary, and looking up (or
 * interpolating) the HEALPix pixel. Off-sky pixels and UNSEEN/BAD_DATA values are NaN.
 */
class CartaHealpixImage : public casacore::ImageInterface<float> {
public:
    CartaHealpixImage(const std::string& filename, unsigned int hdu, const HealpixProjectionOptions& options);
    CartaHealpixImage(const CartaHealpixImage& other);
    ~CartaHealpixImage() override = default;

    // Metadata helpers (no image data read)
    static bool IsHealpixHdu(const std::string& filename, unsigned int hdu, HealpixMapInfo* info = nullptr);
    static bool IsHealpixFile(const std::string& filename, unsigned int* first_hdu = nullptr);
    static bool ReadMapInfo(const std::string& filename, unsigned int hdu, HealpixMapInfo& info, std::string& error);

    // implement casacore ImageInterface
    casacore::String imageType() const override;
    casacore::String name(bool stripPath = false) const override;
    casacore::IPosition shape() const override;
    casacore::Bool ok() const override;
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

    const HealpixMapInfo& MapInfo() const {
        return _info;
    }
    const HealpixProjectionOptions& Options() const {
        return _options;
    }
    std::vector<std::string> ColumnNames() const {
        return _column_names;
    }
    /// Number of maps along the third axis (table columns selected, or image-format maps)
    size_t NumMaps() const;

private:
    void ReadMap();
    void ReadImageFormatMaps();
    void SetUpCoordinateSystem();
    bool AddMapAxisCoordinate(casacore::CoordinateSystem& coord_sys);
    void SetUpFrameRotation();
    void Lookup(const double* lon, const double* lat, const std::vector<bool>& valid, size_t npoints, size_t column, float* out) const;

    std::string _filename;
    unsigned int _hdu;
    HealpixProjectionOptions _options;
    HealpixMapInfo _info;

    std::vector<std::string> _column_names;    // selected columns, in axis order
    std::shared_ptr<HealpixMapSource> _source; // map values (NaN where unseen)
    std::shared_ptr<T_Healpix_Base<int64_t>> _healpix;

    casacore::IPosition _shape;
    casacore::DirectionCoordinate _direction_coord;
    casacore::MDirection::Types _grid_frame; // frame of the direction coordinate
    casacore::MDirection::Types _map_frame;  // frame of the HEALPix map
    bool _needs_rotation;
    casacore::Matrix<double> _rotation; // 3x3, grid frame -> map frame
};

} // namespace carta

#endif // CARTA_SRC_IMAGEDATA_CARTAHEALPIXIMAGE_H_
