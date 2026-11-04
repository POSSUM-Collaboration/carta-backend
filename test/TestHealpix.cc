/* This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-backend
   Copyright 2018- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
   Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
   SPDX-License-Identifier: GPL-3.0-or-later
*/

#include <cmath>

#include <fitsio.h>
#include <gtest/gtest.h>

#include <casacore/casa/Quanta/MVDirection.h>
#include <casacore/measures/Measures/MCDirection.h>
#include <casacore/measures/Measures/MeasConvert.h>

#include <healpix_base.h>
#include <pointing.h>

#include "CommonTestUtilities.h"
#include "FileList/FileExtInfoLoader.h"
#include "FileList/FitsHduList.h"
#include "Frame/Frame.h"
#include "ImageData/CartaHealpixImage.h"
#include "ImageData/FileLoader.h"
#include "ImageData/HealpixLoader.h"

using namespace carta;

namespace {

const double DEG = 180.0 / M_PI;

struct HealpixColumn {
    std::string name;
    std::string unit;
    std::vector<float> values;
};

// Write a HEALPix FITS file (empty primary HDU + binary table extension) with cfitsio
void WriteHealpixFile(const fs::path& path, long nside, bool nested, const std::string& coordsys, const std::vector<HealpixColumn>& columns,
    const std::vector<int64_t>* explicit_pixels = nullptr) {
    fs::remove(path);
    fitsfile* fptr(nullptr);
    int status(0);
    fits_create_file(&fptr, path.string().c_str(), &status);
    ASSERT_EQ(status, 0);

    // Empty primary
    fits_create_img(fptr, 8, 0, nullptr, &status);
    ASSERT_EQ(status, 0);

    // Binary table
    int ncols = columns.size() + (explicit_pixels ? 1 : 0);
    std::vector<std::string> ttype, tform, tunit;
    if (explicit_pixels) {
        ttype.push_back("PIXEL");
        tform.push_back("1K");
        tunit.push_back("");
    }
    for (auto& column : columns) {
        ttype.push_back(column.name);
        tform.push_back("1E");
        tunit.push_back(column.unit);
    }
    std::vector<char*> ttype_ptr, tform_ptr, tunit_ptr;
    for (int i = 0; i < ncols; ++i) {
        ttype_ptr.push_back(const_cast<char*>(ttype[i].c_str()));
        tform_ptr.push_back(const_cast<char*>(tform[i].c_str()));
        tunit_ptr.push_back(const_cast<char*>(tunit[i].c_str()));
    }
    long nrows = explicit_pixels ? explicit_pixels->size() : columns[0].values.size();
    char extname[] = "xtension";
    fits_create_tbl(fptr, BINARY_TBL, nrows, ncols, ttype_ptr.data(), tform_ptr.data(), tunit_ptr.data(), extname, &status);
    ASSERT_EQ(status, 0);

    char pixtype[] = "HEALPIX";
    fits_write_key(fptr, TSTRING, "PIXTYPE", pixtype, "HEALPIX pixelisation", &status);
    std::string ordering_str = nested ? "NESTED" : "RING";
    fits_write_key(fptr, TSTRING, "ORDERING", const_cast<char*>(ordering_str.c_str()), "Pixel ordering scheme", &status);
    fits_write_key(fptr, TLONG, "NSIDE", &nside, "Resolution parameter", &status);
    if (!coordsys.empty()) {
        fits_write_key(fptr, TSTRING, "COORDSYS", const_cast<char*>(coordsys.c_str()), "Coordinate system", &status);
    }
    std::string indxschm = explicit_pixels ? "EXPLICIT" : "IMPLICIT";
    fits_write_key(fptr, TSTRING, "INDXSCHM", const_cast<char*>(indxschm.c_str()), "Indexing scheme", &status);
    std::string object = explicit_pixels ? "PARTIAL" : "FULLSKY";
    fits_write_key(fptr, TSTRING, "OBJECT", const_cast<char*>(object.c_str()), "Sky coverage", &status);
    ASSERT_EQ(status, 0);

    int colnum = 1;
    if (explicit_pixels) {
        fits_write_col(fptr, TLONGLONG, colnum++, 1, 1, nrows, const_cast<int64_t*>(explicit_pixels->data()), &status);
    }
    for (auto& column : columns) {
        fits_write_col(fptr, TFLOAT, colnum++, 1, 1, nrows, const_cast<float*>(column.values.data()), &status);
    }
    ASSERT_EQ(status, 0);
    fits_close_file(fptr, &status);
    ASSERT_EQ(status, 0);
}

// Columns holding the longitude and latitude (deg) of each pixel centre
std::vector<HealpixColumn> LonLatColumns(long nside, bool nested) {
    T_Healpix_Base<int64_t> base(nside, nested ? NEST : RING, SET_NSIDE);
    long npix = 12 * nside * nside;
    HealpixColumn lon{"LON", "deg", std::vector<float>(npix)};
    HealpixColumn lat{"LAT", "deg", std::vector<float>(npix)};
    for (long pix = 0; pix < npix; ++pix) {
        pointing ptg = base.pix2ang(pix);
        lon.values[pix] = ptg.phi * DEG;
        lat.values[pix] = (M_PI / 2.0 - ptg.theta) * DEG;
    }
    return {lon, lat};
}

// Write an image-format HEALPix cube: primary HDU [npix, nmaps] (one full-sky map per row), map k = LON + 1000 k,
// with a FREQ WCS on the map axis numbered as a cube's third axis (as GMIMS does), a beam and a unit
void WriteHealpixImageCube(const fs::path& path, long nside, int nmaps) {
    fitsfile* fptr(nullptr);
    int status(0);
    fs::remove(path);
    fits_create_file(&fptr, path.string().c_str(), &status);
    long npix = 12 * nside * nside;
    long naxes[2] = {npix, nmaps};
    fits_create_img(fptr, FLOAT_IMG, 2, naxes, &status);
    char pixtype[] = "HEALPIX", ordering[] = "RING", coordsys[] = "C", indxschm[] = "IMPLICIT", object[] = "FULLSKY";
    char ctype3[] = "FREQ", cunit3[] = "Hz", bunit[] = "K";
    fits_write_key(fptr, TSTRING, "PIXTYPE", pixtype, nullptr, &status);
    fits_write_key(fptr, TSTRING, "ORDERING", ordering, nullptr, &status);
    fits_write_key(fptr, TLONG, "NSIDE", &nside, nullptr, &status);
    fits_write_key(fptr, TSTRING, "COORDSYS", coordsys, nullptr, &status);
    fits_write_key(fptr, TSTRING, "INDXSCHM", indxschm, nullptr, &status);
    fits_write_key(fptr, TSTRING, "OBJECT", object, nullptr, &status);
    double crval3(286.25e6), cdelt3(0.5e6), crpix3(1.0), bmaj(1.35), bmin(1.2), bpa(10.0), restfreq(300.0e6);
    fits_write_key(fptr, TSTRING, "CTYPE3", ctype3, nullptr, &status);
    fits_write_key(fptr, TDOUBLE, "CRVAL3", &crval3, nullptr, &status);
    fits_write_key(fptr, TDOUBLE, "CDELT3", &cdelt3, nullptr, &status);
    fits_write_key(fptr, TDOUBLE, "CRPIX3", &crpix3, nullptr, &status);
    fits_write_key(fptr, TSTRING, "CUNIT3", cunit3, nullptr, &status);
    fits_write_key(fptr, TDOUBLE, "BMAJ", &bmaj, nullptr, &status);
    fits_write_key(fptr, TDOUBLE, "BMIN", &bmin, nullptr, &status);
    fits_write_key(fptr, TDOUBLE, "BPA", &bpa, nullptr, &status);
    fits_write_key(fptr, TDOUBLE, "RESTFREQ", &restfreq, nullptr, &status);
    fits_write_key(fptr, TSTRING, "BUNIT", bunit, nullptr, &status);

    auto lon = LonLatColumns(nside, false)[0].values;
    std::vector<float> values(npix * nmaps);
    for (int k = 0; k < nmaps; ++k) {
        for (long pix = 0; pix < npix; ++pix) {
            values[k * npix + pix] = lon[pix] + 1000.0f * k;
        }
    }
    values[2] = -1.6375e30f; // an UNSEEN pixel in map 0
    fits_write_img(fptr, TFLOAT, 1, npix * nmaps, values.data(), &status);
    fits_close_file(fptr, &status);
    ASSERT_EQ(status, 0) << "writing " << path;
}

double LonDiff(double a, double b) {
    double d = std::fmod(std::abs(a - b), 360.0);
    return std::min(d, 360.0 - d);
}

// Great-circle separation in degrees between two (lon, lat) positions in degrees
double AngularSeparation(double lon1, double lat1, double lon2, double lat2) {
    double a1 = lon1 / DEG, b1 = lat1 / DEG, a2 = lon2 / DEG, b2 = lat2 / DEG;
    double cos_sep = std::sin(b1) * std::sin(b2) + std::cos(b1) * std::cos(b2) * std::cos(a1 - a2);
    return std::acos(std::max(-1.0, std::min(1.0, cos_sep))) * DEG;
}

fs::path GeneratedPath(const std::string& name) {
    auto dir = fs::temp_directory_path() / "carta_healpix_tests";
    fs::create_directories(dir);
    return dir / name;
}

} // namespace

class HealpixTest : public ::testing::Test {
public:
    static void SetUpTestSuite() {
        WriteHealpixFile(GeneratedPath("lonlat_ring_c.fits"), 32, false, "C", LonLatColumns(32, false));
        WriteHealpixFile(GeneratedPath("lonlat_nest_g.fits"), 16, true, "G", LonLatColumns(16, true));

        // I/Q/U map, nested, galactic
        long nside(16), npix(12 * nside * nside);
        std::vector<HealpixColumn> iqu = {{"I_STOKES", "K_CMB", std::vector<float>(npix, 1.0f)},
            {"Q_STOKES", "K_CMB", std::vector<float>(npix, 2.0f)}, {"U_STOKES", "K_CMB", std::vector<float>(npix, 3.0f)}};
        WriteHealpixFile(GeneratedPath("iqu_nest_g.fits"), nside, true, "G", iqu);

        // Partial map: only the northern hemisphere pixels of a ring map, value = pixel index
        T_Healpix_Base<int64_t> base(32, RING, SET_NSIDE);
        std::vector<int64_t> pixels;
        std::vector<float> values;
        for (long pix = 0; pix < 12 * 32 * 32; ++pix) {
            if (base.pix2ang(pix).theta < M_PI / 2.0) {
                pixels.push_back(pix);
                values.push_back(pix);
            }
        }
        WriteHealpixFile(GeneratedPath("partial_ring_c.fits"), 32, false, "C", {{"SIGNAL", "", values}}, &pixels);

        // Image-format cubes: plain (memory-mapped) and gzip-compressed (read through cfitsio)
        WriteHealpixImageCube(GeneratedPath("cube_image_format.fits"), 16, 3);
        WriteHealpixImageCube(GeneratedPath("cube_image_format.fits.gz"), 16, 3);
    }
};

TEST_F(HealpixTest, DetectsHealpixFile) {
    auto path = GeneratedPath("lonlat_ring_c.fits").string();
    unsigned int hdu(0);
    EXPECT_TRUE(CartaHealpixImage::IsHealpixFile(path, &hdu));
    EXPECT_EQ(hdu, 1);
    EXPECT_FALSE(CartaHealpixImage::IsHealpixFile((FitsImages() / "nonexistent.fits").string()));

    HealpixMapInfo info;
    std::string error;
    ASSERT_TRUE(CartaHealpixImage::ReadMapInfo(path, 1, info, error)) << error;
    EXPECT_EQ(info.nside, 32);
    EXPECT_FALSE(info.nested);
    EXPECT_FALSE(info.explicit_index);
    EXPECT_EQ(info.coordsys, "C");
    ASSERT_EQ(info.columns.size(), 2);
    EXPECT_EQ(info.columns[0], "LON");
    EXPECT_EQ(info.column_units[1], "DEG");
    EXPECT_NEAR(info.resolution_deg(), 1.8323, 1e-3);
}

TEST_F(HealpixTest, HduListAndLoaderDispatch) {
    auto path = GeneratedPath("lonlat_ring_c.fits").string();
    std::vector<std::string> hdu_list;
    std::string error;
    FitsHduList(path).GetHduList(hdu_list, error);
    ASSERT_EQ(hdu_list.size(), 1) << error;
    EXPECT_EQ(hdu_list[0], "1: HEALPix NSIDE=32 RING");

    auto loader = FileLoader::GetLoader(path);
    ASSERT_TRUE(loader);
    EXPECT_TRUE(std::dynamic_pointer_cast<HealpixLoader>(loader) != nullptr);
}

TEST_F(HealpixTest, DefaultMollweideGrid) {
    auto path = GeneratedPath("lonlat_ring_c.fits").string();
    auto loader = FileLoader::GetLoader(path);
    loader->OpenFile("1");
    auto image = loader->GetImage();
    ASSERT_TRUE(image);
    EXPECT_EQ(image->imageType(), "CartaHealpixImage");

    // Default grid: all-sky at the HEALPix resolution, columns on a linear axis
    auto shape = image->shape();
    ASSERT_EQ(shape.size(), 3);
    EXPECT_EQ(shape(0), static_cast<int>(std::ceil(360.0 / 1.8323)));
    EXPECT_EQ(shape(1), static_cast<int>(std::ceil(180.0 / 1.8323)));
    EXPECT_EQ(shape(2), 2);
    EXPECT_TRUE(image->coordinates().hasDirectionCoordinate());
    EXPECT_FALSE(image->coordinates().hasPolarizationAxis());
    EXPECT_EQ(image->coordinates().directionCoordinate().projection().name(), "MOL");

    std::string message;
    ASSERT_TRUE(loader->FindCoordinateAxes(message)) << message;
    auto axes = loader->GetAxes();
    EXPECT_EQ(axes.x, 0);
    EXPECT_EQ(axes.y, 1);
    EXPECT_EQ(axes.z, 2);
    EXPECT_EQ(axes.stokes, -1);

    // Values: the LON and LAT columns give the centre of the HEALPix pixel sampled at each grid pixel,
    // which must lie within a HEALPix pixel size of the grid pixel's own sky position
    casacore::Array<float> data;
    casacore::IPosition start(3, 0, 0, 0), end(3, shape(0) - 1, shape(1) - 1, 1);
    image->doGetSlice(data, casacore::Slicer(start, end, casacore::Slicer::endIsLast));
    const auto& dir = image->coordinates().directionCoordinate();
    int on_sky(0), off_sky(0);
    for (int y = 0; y < shape(1); y += 7) {
        for (int x = 0; x < shape(0); x += 7) {
            float lon_value = data(casacore::IPosition(3, x, y, 0));
            float lat_value = data(casacore::IPosition(3, x, y, 1));
            casacore::Vector<double> pixel(2), world(2);
            pixel(0) = x;
            pixel(1) = y;
            if (dir.toWorld(world, pixel)) {
                EXPECT_TRUE(std::isfinite(lon_value) && std::isfinite(lat_value)) << "pixel " << x << "," << y;
                EXPECT_LT(AngularSeparation(lon_value, lat_value, world(0) * DEG, world(1) * DEG), 1.5 * 1.8323)
                    << "pixel " << x << "," << y;
                ++on_sky;
            } else {
                EXPECT_TRUE(std::isnan(lon_value)) << "pixel " << x << "," << y;
                ++off_sky;
            }
        }
    }
    EXPECT_GT(on_sky, 100);
    EXPECT_GT(off_sky, 10); // Mollweide corners are off the sky
}

TEST_F(HealpixTest, StokesColumnsAndFrame) {
    auto path = GeneratedPath("iqu_nest_g.fits").string();
    auto loader = FileLoader::GetLoader(path);
    loader->OpenFile("");
    auto image = loader->GetImage();
    ASSERT_TRUE(image);

    auto shape = image->shape();
    ASSERT_EQ(shape.size(), 3);
    EXPECT_EQ(shape(2), 3);
    ASSERT_TRUE(image->coordinates().hasPolarizationAxis());
    auto stokes = image->coordinates().stokesCoordinate().stokesStrings();
    ASSERT_EQ(stokes.size(), 3);
    EXPECT_EQ(stokes(0), "I");
    EXPECT_EQ(stokes(1), "Q");
    EXPECT_EQ(stokes(2), "U");
    EXPECT_EQ(image->coordinates().directionCoordinate().directionType(), casacore::MDirection::GALACTIC);

    std::string message;
    ASSERT_TRUE(loader->FindCoordinateAxes(message)) << message;
    EXPECT_EQ(loader->GetAxes().stokes, 2);
    EXPECT_EQ(loader->GetDims().num_stokes, 3);

    // Each Stokes plane is constant
    casacore::Array<float> data;
    int cx = shape(0) / 2, cy = shape(1) / 2;
    casacore::IPosition start(3, cx, cy, 0), end(3, cx, cy, 2);
    image->doGetSlice(data, casacore::Slicer(start, end, casacore::Slicer::endIsLast));
    EXPECT_FLOAT_EQ(data(casacore::IPosition(3, 0, 0, 0)), 1.0f);
    EXPECT_FLOAT_EQ(data(casacore::IPosition(3, 0, 0, 1)), 2.0f);
    EXPECT_FLOAT_EQ(data(casacore::IPosition(3, 0, 0, 2)), 3.0f);
}

TEST_F(HealpixTest, ProjectionOptions) {
    auto path = GeneratedPath("lonlat_nest_g.fits").string();
    auto loader = std::dynamic_pointer_cast<HealpixLoader>(FileLoader::GetLoader(path));
    ASSERT_TRUE(loader);

    HealpixProjectionOptions options;
    options.projection = "TAN";
    options.frame = "C"; // reproject the galactic map onto an equatorial tangent plane
    options.center_lon_deg = 266.4;
    options.center_lat_deg = -28.9; // near the galactic centre
    options.pixel_scale_deg = 0.5;
    options.width = 40;
    options.height = 30;
    options.columns = {"LAT"};
    loader->SetProjectionOptions(options);
    loader->OpenFile("");
    auto image = loader->GetImage();
    ASSERT_TRUE(image);

    auto shape = image->shape();
    ASSERT_EQ(shape.size(), 2);
    EXPECT_EQ(shape(0), 40);
    EXPECT_EQ(shape(1), 30);
    const auto& dir = image->coordinates().directionCoordinate();
    EXPECT_EQ(dir.projection().name(), "TAN");
    EXPECT_EQ(dir.directionType(), casacore::MDirection::ICRS);

    // The LAT column holds galactic latitude: convert each grid pixel ICRS -> Galactic and compare
    casacore::Array<float> data;
    casacore::IPosition start(2, 0, 0), end(2, 39, 29);
    image->doGetSlice(data, casacore::Slicer(start, end, casacore::Slicer::endIsLast));
    casacore::MDirection::Ref icrs_ref(casacore::MDirection::ICRS), galactic_ref(casacore::MDirection::GALACTIC);
    casacore::MDirection::Convert to_galactic(icrs_ref, galactic_ref);
    const double resolution = 3.66; // deg at NSIDE=16
    for (int y = 0; y < 30; y += 5) {
        for (int x = 0; x < 40; x += 5) {
            casacore::Vector<double> pixel(2), world(2);
            pixel(0) = x;
            pixel(1) = y;
            ASSERT_TRUE(dir.toWorld(world, pixel));
            casacore::MDirection icrs(casacore::Quantity(world(0), "rad"), casacore::Quantity(world(1), "rad"), casacore::MDirection::ICRS);
            auto galactic = to_galactic(icrs).getValue().get();
            double glat = galactic(1) * DEG;
            float value = data(casacore::IPosition(2, x, y));
            ASSERT_TRUE(std::isfinite(value));
            EXPECT_NEAR(value, glat, 2.0 * resolution) << "pixel " << x << "," << y;
        }
    }
}

TEST_F(HealpixTest, MatchToImageGrid) {
    auto path = GeneratedPath("lonlat_ring_c.fits").string();
    auto loader = std::dynamic_pointer_cast<HealpixLoader>(FileLoader::GetLoader(path));
    ASSERT_TRUE(loader);

    // Target: a 4D [RA, DEC, FREQ, STOKES]-like system with a SIN projection, lat axis first
    casacore::Matrix<double> xform(2, 2, 0.0);
    xform.diagonal() = 1.0;
    casacore::DirectionCoordinate target_dir(casacore::MDirection::J2000, casacore::Projection(casacore::Projection::SIN), 80.0 / DEG,
        -40.0 / DEG, -0.25 / DEG, 0.25 / DEG, xform, 30.0, 20.0);
    auto match = std::make_shared<casacore::CoordinateSystem>();
    match->addCoordinate(target_dir);
    casacore::Vector<casacore::Int> stokes_types(1, casacore::Stokes::I);
    match->addCoordinate(casacore::StokesCoordinate(stokes_types));
    casacore::IPosition match_shape(3, 60, 40, 1);

    HealpixProjectionOptions options;
    options.match_coord_sys = match;
    options.match_shape = match_shape;
    options.columns = {"LON"};
    options.bilinear = true;
    loader->SetProjectionOptions(options);
    loader->OpenFile("1");
    auto image = loader->GetImage();
    ASSERT_TRUE(image);

    auto shape = image->shape();
    ASSERT_EQ(shape.size(), 2);
    EXPECT_EQ(shape(0), 60);
    EXPECT_EQ(shape(1), 40);
    const auto& dir = image->coordinates().directionCoordinate();
    EXPECT_EQ(dir.projection().name(), "SIN");
    casacore::Vector<double> refval = dir.referenceValue();
    EXPECT_NEAR(refval(0) * DEG, 80.0, 1e-9);
    EXPECT_NEAR(refval(1) * DEG, -40.0, 1e-9);

    casacore::Array<float> data;
    casacore::IPosition start(2, 0, 0), end(2, 59, 39);
    image->doGetSlice(data, casacore::Slicer(start, end, casacore::Slicer::endIsLast));
    for (int y = 0; y < 40; y += 8) {
        for (int x = 0; x < 60; x += 8) {
            casacore::Vector<double> pixel(2), world(2);
            pixel(0) = x;
            pixel(1) = y;
            ASSERT_TRUE(dir.toWorld(world, pixel));
            float value = data(casacore::IPosition(2, x, y));
            ASSERT_TRUE(std::isfinite(value));
            EXPECT_LT(LonDiff(value, world(0) * DEG), 2.0 * 1.8323) << "pixel " << x << "," << y;
        }
    }

    // Frame opens it like any other image
    std::unique_ptr<Frame> frame(new Frame(0, loader, "1"));
    EXPECT_TRUE(frame->IsValid()) << frame->GetErrorMessage();
    EXPECT_EQ(frame->Width(), 60);
    EXPECT_EQ(frame->Height(), 40);
    EXPECT_EQ(frame->Depth(), 1);
}

TEST_F(HealpixTest, PartialMapHasNaNOutsideCoverage) {
    auto path = GeneratedPath("partial_ring_c.fits").string();
    auto loader = std::dynamic_pointer_cast<HealpixLoader>(FileLoader::GetLoader(path));
    ASSERT_TRUE(loader);

    HealpixProjectionOptions options;
    options.projection = "CAR";
    options.pixel_scale_deg = 2.0;
    options.width = 180;
    options.height = 90;
    loader->SetProjectionOptions(options);
    loader->OpenFile("");
    auto image = loader->GetImage();
    ASSERT_TRUE(image);
    EXPECT_EQ(loader->GetImage()->shape()(0), 180);

    casacore::Array<float> data;
    casacore::IPosition start(2, 0, 0), end(2, 179, 89);
    image->doGetSlice(data, casacore::Slicer(start, end, casacore::Slicer::endIsLast));
    const auto& dir = image->coordinates().directionCoordinate();
    int seen(0), unseen(0);
    for (int y = 0; y < 90; y += 3) {
        casacore::Vector<double> pixel(2), world(2);
        pixel(0) = 90;
        pixel(1) = y;
        ASSERT_TRUE(dir.toWorld(world, pixel));
        float value = data(casacore::IPosition(2, 90, y));
        if (world(1) > 0.03) { // comfortably in the northern hemisphere
            EXPECT_TRUE(std::isfinite(value)) << "row " << y;
            ++seen;
        } else if (world(1) < -0.03) {
            EXPECT_TRUE(std::isnan(value)) << "row " << y;
            ++unseen;
        }
    }
    EXPECT_GT(seen, 10);
    EXPECT_GT(unseen, 10);
}

TEST_F(HealpixTest, FileInfoEntries) {
    auto path = GeneratedPath("iqu_nest_g.fits").string();
    auto loader = FileLoader::GetLoader(path);
    FileExtInfoLoader ext_info_loader(loader);
    std::map<std::string, CARTA::FileInfoExtended> info_map;
    std::string message;
    ASSERT_TRUE(ext_info_loader.FillFitsFileInfoMap(info_map, path, message)) << message;
    ASSERT_EQ(info_map.size(), 1);
    const auto& info = info_map.begin()->second;
    EXPECT_EQ(info.dimensions(), 3);
    EXPECT_EQ(info.stokes(), 3);

    std::map<std::string, std::string> computed;
    for (int i = 0; i < info.computed_entries_size(); ++i) {
        computed[info.computed_entries(i).name()] = info.computed_entries(i).value();
    }
    EXPECT_EQ(computed["HEALPix NSIDE"], "16");
    EXPECT_EQ(computed["HEALPix ordering"], "NESTED");
    EXPECT_EQ(computed["HEALPix frame"], "G");
    EXPECT_EQ(computed["HEALPix columns"], "I_STOKES [K_CMB], Q_STOKES [K_CMB], U_STOKES [K_CMB]");
    EXPECT_NE(computed["HEALPix reprojection"].find("MOL"), std::string::npos);
    EXPECT_EQ(computed["Projection"], "MOL");
    EXPECT_EQ(computed["Number of polarizations"], "3");
    // The synthesised header must carry the WCS for the frontend
    bool has_ctype1(false);
    for (int i = 0; i < info.header_entries_size(); ++i) {
        if (info.header_entries(i).name() == "CTYPE1") {
            has_ctype1 = true;
            EXPECT_NE(info.header_entries(i).value().find("GLON-MOL"), std::string::npos);
        }
    }
    EXPECT_TRUE(has_ctype1);
}

TEST_F(HealpixTest, LargeSectionIsReadInChunksConsistently) {
    // A section larger than the chunk size (2^18 points), with a stride and an offset, must equal the same
    // pixels read in small pieces: chunk boundaries fall at arbitrary (x, y) positions
    auto path = GeneratedPath("lonlat_ring_c.fits").string();
    auto loader = std::dynamic_pointer_cast<HealpixLoader>(FileLoader::GetLoader(path));
    ASSERT_TRUE(loader);
    HealpixProjectionOptions options;
    options.projection = "MOL";
    options.pixel_scale_deg = 0.25;
    options.width = 1440;
    options.height = 720;
    options.columns = {"LON"};
    loader->SetProjectionOptions(options);
    loader->OpenFile("");
    auto image = loader->GetImage();
    ASSERT_TRUE(image);

    casacore::Array<float> whole;
    casacore::IPosition start(2, 3, 7), end(2, 1437, 715), stride(2, 1, 1);
    image->doGetSlice(whole, casacore::Slicer(start, end, stride, casacore::Slicer::endIsLast));
    ASSERT_EQ(whole.shape(), casacore::IPosition(2, 1435, 709));
    ASSERT_GT(whole.shape().product(), 1 << 18);

    size_t finite(0);
    for (int y = 0; y < 709; y += 97) {
        casacore::Array<float> row;
        casacore::IPosition row_start(2, 3, 7 + y), row_end(2, 1437, 7 + y);
        image->doGetSlice(row, casacore::Slicer(row_start, row_end, casacore::Slicer::endIsLast));
        for (int x = 0; x < 1435; x += 13) {
            float a = whole(casacore::IPosition(2, x, y));
            float b = row(casacore::IPosition(2, x, 0));
            if (std::isfinite(a) || std::isfinite(b)) {
                ASSERT_TRUE(std::isfinite(a) && std::isfinite(b)) << "pixel " << x << "," << y;
                ASSERT_FLOAT_EQ(a, b) << "pixel " << x << "," << y;
                ++finite;
            }
        }
    }
    EXPECT_GT(finite, 500);

    // A strided read picks the matching pixels of the full read
    casacore::Array<float> strided;
    casacore::IPosition stride3(2, 3, 2);
    image->doGetSlice(strided, casacore::Slicer(start, end, stride3, casacore::Slicer::endIsLast));
    for (int y = 0; y < static_cast<int>(strided.shape()(1)); y += 50) {
        for (int x = 0; x < static_cast<int>(strided.shape()(0)); x += 50) {
            float a = whole(casacore::IPosition(2, 3 * x, 2 * y));
            float b = strided(casacore::IPosition(2, x, y));
            if (std::isfinite(a) || std::isfinite(b)) {
                ASSERT_FLOAT_EQ(a, b) << "strided pixel " << x << "," << y;
            }
        }
    }
}

TEST_F(HealpixTest, ImageFormatCubeOpensAsSpectralCube) {
    for (const std::string name : {"cube_image_format.fits", "cube_image_format.fits.gz"}) {
        auto path = GeneratedPath(name).string();
        unsigned int first_hdu(99);
        ASSERT_TRUE(CartaHealpixImage::IsHealpixFile(path, &first_hdu)) << name;
        EXPECT_EQ(first_hdu, 0);

        HealpixMapInfo info;
        std::string error;
        ASSERT_TRUE(CartaHealpixImage::ReadMapInfo(path, 0, info, error)) << error;
        EXPECT_TRUE(info.image_format);
        EXPECT_EQ(info.num_maps, 3);
        EXPECT_EQ(info.nside, 16);
        EXPECT_EQ(info.map_axis_ctype, "FREQ");
        EXPECT_TRUE(info.has_map_axis());
        EXPECT_DOUBLE_EQ(info.map_axis_crval, 286.25e6);
        EXPECT_EQ(info.columns, std::vector<std::string>({"MAP"}));
        EXPECT_EQ(info.column_units, std::vector<std::string>({"K"}));

        if (name == "cube_image_format.fits") { // gzipped files take the CompressedFits header path instead
            std::vector<std::string> hdu_list;
            FitsHduList fits_hdu_list(path);
            fits_hdu_list.GetHduList(hdu_list, error);
            ASSERT_EQ(hdu_list.size(), 1) << error;
            EXPECT_EQ(hdu_list[0], "0: HEALPix NSIDE=16 RING");
        }

        auto loader = std::dynamic_pointer_cast<HealpixLoader>(FileLoader::GetLoader(path));
        ASSERT_TRUE(loader) << name;
        HealpixProjectionOptions options;
        options.pixel_scale_deg = 2.0;
        loader->SetProjectionOptions(options);
        loader->OpenFile("");
        auto image = loader->GetImage();
        ASSERT_TRUE(image);
        auto shape = image->shape();
        ASSERT_EQ(shape.size(), 3) << name;
        EXPECT_EQ(shape(0), 180);
        EXPECT_EQ(shape(1), 90);
        EXPECT_EQ(shape(2), 3);
        EXPECT_EQ(image->units().getName(), "K");
        EXPECT_FALSE(image->imageInfo().restoringBeam().isNull());
        EXPECT_NEAR(image->imageInfo().restoringBeam().getMajor("deg"), 1.35, 1e-9);

        const auto& coord_sys = image->coordinates();
        ASSERT_TRUE(coord_sys.hasSpectralAxis());
        EXPECT_EQ(coord_sys.spectralAxisNumber(), 2);
        double freq(0.0);
        ASSERT_TRUE(coord_sys.spectralCoordinate().toWorld(freq, 2.0));
        EXPECT_NEAR(freq, 287.25e6, 1.0);
        EXPECT_NEAR(coord_sys.spectralCoordinate().restFrequency(), 300.0e6, 1.0);

        // Map k holds LON + 1000 k: read channel 2 and compare with the grid's own longitudes
        casacore::Array<float> data;
        casacore::IPosition start(3, 0, 0, 2), end(3, 179, 89, 2);
        image->doGetSlice(data, casacore::Slicer(start, end, casacore::Slicer::endIsLast));
        const auto& dir = coord_sys.directionCoordinate();
        const double resolution = 3.66; // deg at NSIDE=16
        int checked(0);
        for (int y = 5; y < 90; y += 10) {
            for (int x = 5; x < 180; x += 15) {
                casacore::Vector<double> pixel(2), world(2);
                pixel(0) = x;
                pixel(1) = y;
                if (!dir.toWorld(world, pixel)) {
                    continue;
                }
                float value = data(casacore::IPosition(3, x, y, 0));
                ASSERT_TRUE(std::isfinite(value)) << name << " pixel " << x << "," << y;
                // longitude difference scaled to an angle on the sky (pixels are wide in longitude near the poles)
                double lon_error = LonDiff(value - 2000.0f, world(0) * DEG) * std::cos(world(1));
                EXPECT_LT(lon_error, 2.0 * resolution) << name << " pixel " << x << "," << y;
                ++checked;
            }
        }
        EXPECT_GT(checked, 50);

        // A spectrum through one position reads every map; map 0 has an UNSEEN pixel that becomes NaN
        casacore::Array<float> spectrum;
        casacore::IPosition spec_start(3, 90, 45, 0), spec_end(3, 90, 45, 2);
        image->doGetSlice(spectrum, casacore::Slicer(spec_start, spec_end, casacore::Slicer::endIsLast));
        ASSERT_EQ(spectrum.shape(), casacore::IPosition(3, 1, 1, 3));
        float v0 = spectrum(casacore::IPosition(3, 0, 0, 0)), v1 = spectrum(casacore::IPosition(3, 0, 0, 1)),
              v2 = spectrum(casacore::IPosition(3, 0, 0, 2));
        EXPECT_NEAR(v1 - v0, 1000.0f, 1e-3);
        EXPECT_NEAR(v2 - v1, 1000.0f, 1e-3);

        // Frame treats it as a spectral cube
        std::unique_ptr<Frame> frame(new Frame(0, loader, "0"));
        ASSERT_TRUE(frame->IsValid()) << frame->GetErrorMessage();
        EXPECT_EQ(frame->Depth(), 3);
        EXPECT_EQ(frame->Width(), 180);
    }
}

TEST_F(HealpixTest, MatchedGridCanBeEnlargedAroundTheImage) {
    auto path = GeneratedPath("lonlat_ring_c.fits").string();
    auto loader = std::dynamic_pointer_cast<HealpixLoader>(FileLoader::GetLoader(path));
    ASSERT_TRUE(loader);

    casacore::Matrix<double> xform(2, 2, 0.0);
    xform.diagonal() = 1.0;
    casacore::DirectionCoordinate target_dir(casacore::MDirection::J2000, casacore::Projection(casacore::Projection::TAN), 80.0 / DEG,
        -40.0 / DEG, -0.25 / DEG, 0.25 / DEG, xform, 30.0, 20.0);
    auto match = std::make_shared<casacore::CoordinateSystem>();
    match->addCoordinate(target_dir);

    HealpixProjectionOptions options;
    options.match_coord_sys = match;
    options.match_shape = casacore::IPosition(2, 60, 40);
    options.width = 100; // 20 more pixels on each side
    options.height = 60; // 10 more on each side
    options.columns = {"LON"};
    loader->SetProjectionOptions(options);
    loader->OpenFile("1");
    auto image = loader->GetImage();
    ASSERT_TRUE(image);
    EXPECT_EQ(image->shape(), casacore::IPosition(2, 100, 60));

    // Pixel (x + 20, y + 10) of the enlarged grid is pixel (x, y) of the image
    const auto& dir = image->coordinates().directionCoordinate();
    EXPECT_EQ(dir.projection().name(), "TAN");
    for (int y = 0; y < 40; y += 13) {
        for (int x = 0; x < 60; x += 17) {
            casacore::Vector<double> pixel(2), world(2), target_world(2);
            pixel(0) = x + 20;
            pixel(1) = y + 10;
            ASSERT_TRUE(dir.toWorld(world, pixel));
            pixel(0) = x;
            pixel(1) = y;
            ASSERT_TRUE(target_dir.toWorld(target_world, pixel));
            EXPECT_NEAR(world(0), target_world(0), 1e-12) << x << "," << y;
            EXPECT_NEAR(world(1), target_world(1), 1e-12) << x << "," << y;
        }
    }
}

TEST_F(HealpixTest, ExtendedMatchShowsTheMapAroundTheImage) {
    // A small SIN image at 0.25 deg/pixel: extended, the grid reaches the hemisphere edge (1 rad = 229 pixels each way)
    auto path = GeneratedPath("lonlat_ring_c.fits").string();
    auto loader = std::dynamic_pointer_cast<HealpixLoader>(FileLoader::GetLoader(path));
    ASSERT_TRUE(loader);

    casacore::Matrix<double> xform(2, 2, 0.0);
    xform.diagonal() = 1.0;
    casacore::DirectionCoordinate target_dir(casacore::MDirection::J2000, casacore::Projection(casacore::Projection::SIN), 80.0 / DEG,
        -40.0 / DEG, -0.25 / DEG, 0.25 / DEG, xform, 30.0, 20.0);
    auto match = std::make_shared<casacore::CoordinateSystem>();
    match->addCoordinate(target_dir);

    HealpixProjectionOptions options;
    options.match_coord_sys = match;
    options.match_shape = casacore::IPosition(2, 60, 40);
    options.extend_sky = true;
    options.columns = {"LON"};
    loader->SetProjectionOptions(options);
    loader->OpenFile("1");
    auto image = loader->GetImage();
    ASSERT_TRUE(image);
    auto shape = image->shape();
    // 229 pixels of reach each side of the centre plus the centre itself, then rounded to keep the image on integer pixels
    EXPECT_GE(shape(0), 2 * 229 - 2);
    EXPECT_LE(shape(0), 2 * 229 + 2);
    EXPECT_GE(shape(1), 2 * 229 - 2);
    EXPECT_LE(shape(1), 2 * 229 + 2);
    EXPECT_EQ((shape(0) - 60) % 2, 0);
    EXPECT_EQ((shape(1) - 40) % 2, 0);

    // The image's own pixels sit in the middle with the same world coordinates
    const auto& dir = image->coordinates().directionCoordinate();
    const int dx = (shape(0) - 60) / 2, dy = (shape(1) - 40) / 2;
    casacore::Vector<double> pixel(2), world(2), target_world(2);
    pixel(0) = 30 + dx;
    pixel(1) = 20 + dy;
    ASSERT_TRUE(dir.toWorld(world, pixel));
    pixel(0) = 30;
    pixel(1) = 20;
    ASSERT_TRUE(target_dir.toWorld(target_world, pixel));
    EXPECT_NEAR(world(0), target_world(0), 1e-12);
    EXPECT_NEAR(world(1), target_world(1), 1e-12);

    // Tiles come from the loader at every level without a full-plane cache: the frame opens, a mip-4 tile samples
    // block centres, and a full-resolution chunk matches a direct read
    EXPECT_TRUE(loader->UseTileCache());
    EXPECT_TRUE(loader->HasMip(2));
    EXPECT_TRUE(loader->HasMip(8));
    std::unique_ptr<Frame> frame(new Frame(0, loader, "1"));
    ASSERT_TRUE(frame->IsValid()) << frame->GetErrorMessage();
    EXPECT_EQ(frame->Width(), shape(0));

    // A 64x64 region around the image (the grid's corners are off the hemisphere), on a 4-pixel boundary and inside chunk 0
    std::mutex mutex;
    const int x0 = std::min(448, ((static_cast<int>(shape(0)) / 2 - 32) / 4) * 4),
              y0 = std::min(448, ((static_cast<int>(shape(1)) / 2 - 32) / 4) * 4);
    CARTA::ImageBounds bounds;
    bounds.set_x_min(x0);
    bounds.set_x_max(x0 + 64);
    bounds.set_y_min(y0);
    bounds.set_y_max(y0 + 64);
    std::vector<float> mip_data;
    ASSERT_TRUE(loader->GetDownsampledRasterData(mip_data, 0, 0, bounds, 4, mutex));
    ASSERT_EQ(mip_data.size(), 16 * 16);
    std::vector<float> full_data;
    int chunk_width(0), chunk_height(0);
    ASSERT_TRUE(loader->GetChunk(full_data, chunk_width, chunk_height, 0, 0, 0, 0, mutex));
    EXPECT_EQ(chunk_width, std::min(512, static_cast<int>(shape(0))));
    EXPECT_EQ(chunk_height, std::min(512, static_cast<int>(shape(1))));
    casacore::Array<float> direct;
    casacore::IPosition start(2, x0, y0), end(2, x0 + 63, y0 + 63);
    image->doGetSlice(direct, casacore::Slicer(start, end, casacore::Slicer::endIsLast));
    int compared(0);
    for (int j = 0; j < 16; ++j) {
        for (int i = 0; i < 16; ++i) {
            // block centre of a 4x4 block = full-resolution pixel (4i + 1.5, 4j + 1.5): between pixels, so compare with its neighbours
            float m = mip_data[j * 16 + i];
            float a = direct(casacore::IPosition(2, 4 * i + 1, 4 * j + 1)), b = direct(casacore::IPosition(2, 4 * i + 2, 4 * j + 2));
            float c = full_data[(y0 + 4 * j + 1) * chunk_width + x0 + 4 * i + 1];
            if (std::isfinite(m) && std::isfinite(a) && std::isfinite(b)) {
                EXPECT_LE(LonDiff(m, a), 2.0 * 1.8323 + 1e-3) << i << "," << j;
                EXPECT_FLOAT_EQ(a, c);
                ++compared;
            }
        }
    }
    EXPECT_GT(compared, 100);

    // Statistics describe the map itself: LON spans 0..360 over all 12288 pixels
    auto& stats = loader->GetImageStats(0, 0);
    ASSERT_TRUE(stats.valid);
    EXPECT_EQ(stats.basic_stats[CARTA::StatsType::NumPixels], 12 * 32 * 32);
    EXPECT_NEAR(stats.basic_stats[CARTA::StatsType::Min], 0.0, 2.0);
    EXPECT_NEAR(stats.basic_stats[CARTA::StatsType::Max], 360.0, 2.0);
    EXPECT_NEAR(stats.basic_stats[CARTA::StatsType::Mean], 180.0, 1.0);
    EXPECT_GT(stats.histogram_bins.size(), 10);
    long total(0);
    for (auto n : stats.histogram_bins) {
        total += n;
    }
    EXPECT_EQ(total, 12 * 32 * 32);
}
