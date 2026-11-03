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
