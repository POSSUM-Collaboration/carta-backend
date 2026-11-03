/* This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-backend
   Copyright 2018- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
   Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
   SPDX-License-Identifier: GPL-3.0-or-later
*/

// Tests for images with any number of axes: arbitrary axis types, 5D/6D cubes, and pinned extra axes.

#include <cmath>
#include <numeric>

#include <fitsio.h>
#include <gtest/gtest.h>

#include "CommonTestUtilities.h"
#include "FileList/FileExtInfoLoader.h"
#include "Frame/Frame.h"
#include "ImageData/CartaFitsImage.h"
#include "ImageData/FileLoader.h"

using namespace carta;

namespace {

struct Axis {
    std::string ctype;
    std::string cunit;
    long n;
    double crval, cdelt, crpix;
};

// Pixel value encodes the full N-D index: value = sum_k index_k * 16^k (all test axes are shorter than 16),
// which stays exactly representable in a float for up to 6 axes
float EncodedValue(const std::vector<long>& index) {
    double value(0.0);
    double scale(1.0);
    for (auto i : index) {
        value += i * scale;
        scale *= 16.0;
    }
    return static_cast<float>(value);
}

void WriteCube(const fs::path& path, const std::vector<Axis>& axes) {
    fs::remove(path);
    fitsfile* fptr(nullptr);
    int status(0);
    fits_create_file(&fptr, path.string().c_str(), &status);
    ASSERT_EQ(status, 0);

    std::vector<long> naxes;
    for (auto& axis : axes) {
        naxes.push_back(axis.n);
    }
    fits_create_img(fptr, FLOAT_IMG, naxes.size(), naxes.data(), &status);
    ASSERT_EQ(status, 0);

    for (size_t i = 0; i < axes.size(); ++i) {
        std::string n = std::to_string(i + 1);
        fits_write_key(fptr, TSTRING, ("CTYPE" + n).c_str(), const_cast<char*>(axes[i].ctype.c_str()), nullptr, &status);
        if (!axes[i].cunit.empty()) {
            fits_write_key(fptr, TSTRING, ("CUNIT" + n).c_str(), const_cast<char*>(axes[i].cunit.c_str()), nullptr, &status);
        }
        double crval(axes[i].crval), cdelt(axes[i].cdelt), crpix(axes[i].crpix);
        fits_write_key(fptr, TDOUBLE, ("CRVAL" + n).c_str(), &crval, nullptr, &status);
        fits_write_key(fptr, TDOUBLE, ("CDELT" + n).c_str(), &cdelt, nullptr, &status);
        fits_write_key(fptr, TDOUBLE, ("CRPIX" + n).c_str(), &crpix, nullptr, &status);
    }
    char bunit[] = "Jy/beam";
    fits_write_key(fptr, TSTRING, "BUNIT", bunit, nullptr, &status);
    ASSERT_EQ(status, 0);

    long npix = std::accumulate(naxes.begin(), naxes.end(), 1L, std::multiplies<long>());
    std::vector<float> data(npix);
    std::vector<long> index(axes.size(), 0);
    for (long p = 0; p < npix; ++p) {
        data[p] = EncodedValue(index);
        // increment N-D index, first axis fastest
        for (size_t k = 0; k < index.size(); ++k) {
            if (++index[k] < naxes[k]) {
                break;
            }
            index[k] = 0;
        }
    }
    long first(1);
    fits_write_img(fptr, TFLOAT, first, npix, data.data(), &status);
    ASSERT_EQ(status, 0);
    fits_close_file(fptr, &status);
    ASSERT_EQ(status, 0);
}

fs::path GeneratedPath(const std::string& name) {
    auto dir = fs::temp_directory_path() / "carta_nd_tests";
    fs::create_directories(dir);
    return dir / name;
}

const Axis RA{"RA---SIN", "deg", 12, 180.0, -0.01, 6};
const Axis DEC{"DEC--SIN", "deg", 10, -45.0, 0.01, 5};
const Axis FREQ{"FREQ", "Hz", 6, 1.4e9, 1e6, 1};
const Axis STOKES{"STOKES", "", 3, 1, 1, 1};
const Axis FDEP{"FDEP", "rad/m2", 8, -100.0, 10.0, 1};
const Axis TIME{"TIME", "s", 5, 0.0, 1.0, 1};
const Axis BEAM{"BEAM", "", 4, 1, 1, 1};

} // namespace

class NdAxesTest : public ::testing::Test {
public:
    static void SetUpTestSuite() {
        WriteCube(GeneratedPath("ra_dec_freq_stokes_time.fits"), {RA, DEC, FREQ, STOKES, TIME});
        WriteCube(GeneratedPath("ra_dec_fdep_time.fits"), {RA, DEC, FDEP, TIME});
        WriteCube(GeneratedPath("ra_dec_fdep.fits"), {RA, DEC, FDEP});
        WriteCube(GeneratedPath("ra_dec_freq_time_beam_stokes.fits"), {RA, DEC, FREQ, TIME, BEAM, STOKES});
        WriteCube(GeneratedPath("fdep_time_2d.fits"), {FDEP, TIME});
    }

    static void CheckPlane(Frame& frame, const std::vector<long>& expected_index) {
        // Read the current render plane and check that every pixel encodes the expected position on the non-render axes
        int x_axis = frame.XAxis(), y_axis = frame.YAxis();
        auto stokes_slicer = frame.GetImageSlicer(AxisRange(frame.CurrentZ()), frame.CurrentStokes());
        std::vector<float> plane(stokes_slicer.slicer.length().product());
        ASSERT_TRUE(frame.GetSlicerData(stokes_slicer, plane.data()));
        ASSERT_EQ(plane.size(), frame.Width() * frame.Height());

        std::vector<long> index(expected_index);
        for (size_t y = 0; y < frame.Height(); ++y) {
            for (size_t x = 0; x < frame.Width(); ++x) {
                index[x_axis] = x;
                index[y_axis] = y;
                float expected = EncodedValue(index);
                float actual = plane[y * frame.Width() + x];
                ASSERT_FLOAT_EQ(actual, expected) << "pixel " << x << "," << y;
            }
        }
    }
};

TEST_F(NdAxesTest, Cube5dOpens) {
    auto path = GeneratedPath("ra_dec_freq_stokes_time.fits").string();
    auto loader = FileLoader::GetLoader(path);
    ASSERT_TRUE(loader);
    loader->OpenFile("0");
    std::string message;
    ASSERT_TRUE(loader->FindCoordinateAxes(message)) << message;

    auto axes = loader->GetAxes();
    EXPECT_EQ(axes.x, 0);
    EXPECT_EQ(axes.y, 1);
    EXPECT_EQ(axes.z, 2);
    EXPECT_EQ(axes.spectral, 2);
    EXPECT_EQ(axes.stokes, 3);
    ASSERT_EQ(axes.extra, std::vector<int>({4}));

    auto dims = loader->GetDims();
    EXPECT_EQ(dims.width, 12);
    EXPECT_EQ(dims.height, 10);
    EXPECT_EQ(dims.depth, 6);
    EXPECT_EQ(dims.num_stokes, 3);

    Frame frame(0, loader, "0");
    ASSERT_TRUE(frame.IsValid()) << frame.GetErrorMessage();
    EXPECT_EQ(frame.ExtraAxisIndices(), std::vector<int>({0}));
    CheckPlane(frame, {0, 0, 0, 0, 0});

    // Move along z, stokes and the extra TIME axis
    ASSERT_TRUE(frame.SetImageChannels(4, 2, message)) << message;
    CheckPlane(frame, {0, 0, 4, 2, 0});
    ASSERT_TRUE(frame.SetExtraAxisIndices({3}, message)) << message;
    CheckPlane(frame, {0, 0, 4, 2, 3});

    // Invalid indices are rejected
    EXPECT_FALSE(frame.SetExtraAxisIndices({5}, message));
    EXPECT_FALSE(message.empty());
    EXPECT_FALSE(frame.SetExtraAxisIndices({0, 0}, message));
    EXPECT_EQ(frame.ExtraAxisIndices(), std::vector<int>({3}));
}

TEST_F(NdAxesTest, NonStokesFourthAxisIsNotGuessedAsStokes) {
    auto path = GeneratedPath("ra_dec_fdep_time.fits").string();
    auto loader = FileLoader::GetLoader(path);
    loader->OpenFile("0");
    std::string message;
    ASSERT_TRUE(loader->FindCoordinateAxes(message)) << message;

    auto axes = loader->GetAxes();
    EXPECT_EQ(axes.z, 2);
    EXPECT_EQ(axes.spectral, 2); // no spectral coordinate: z is treated as the channel axis
    EXPECT_EQ(axes.stokes, -1);
    EXPECT_EQ(axes.extra, std::vector<int>({3}));
    EXPECT_EQ(loader->GetDims().num_stokes, 1);
    EXPECT_EQ(loader->GetDims().depth, 8);

    Frame frame(0, loader, "0");
    ASSERT_TRUE(frame.IsValid()) << frame.GetErrorMessage();
    ASSERT_TRUE(frame.SetExtraAxisIndices({2}, message)) << message;
    ASSERT_TRUE(frame.SetImageChannels(7, 0, message)) << message;
    CheckPlane(frame, {0, 0, 7, 2});
}

TEST_F(NdAxesTest, FdfCubeUnchanged) {
    auto path = GeneratedPath("ra_dec_fdep.fits").string();
    auto loader = FileLoader::GetLoader(path);
    loader->OpenFile("0");
    std::string message;
    ASSERT_TRUE(loader->FindCoordinateAxes(message)) << message;
    auto axes = loader->GetAxes();
    EXPECT_EQ(axes.z, 2);
    EXPECT_EQ(axes.spectral, 2);
    EXPECT_EQ(axes.stokes, -1);
    EXPECT_TRUE(axes.extra.empty());
}

TEST_F(NdAxesTest, Cube6dWithTrailingStokes) {
    auto path = GeneratedPath("ra_dec_freq_time_beam_stokes.fits").string();
    auto loader = FileLoader::GetLoader(path);
    loader->OpenFile("0");
    std::string message;
    ASSERT_TRUE(loader->FindCoordinateAxes(message)) << message;

    auto axes = loader->GetAxes();
    EXPECT_EQ(axes.z, 2);
    EXPECT_EQ(axes.stokes, 5);
    EXPECT_EQ(axes.extra, std::vector<int>({3, 4}));

    Frame frame(0, loader, "0");
    ASSERT_TRUE(frame.IsValid()) << frame.GetErrorMessage();
    ASSERT_TRUE(frame.SetExtraAxisIndices({4, 1}, message)) << message;
    ASSERT_TRUE(frame.SetImageChannels(5, 1, message)) << message;
    CheckPlane(frame, {0, 0, 5, 4, 1, 1});
}

TEST_F(NdAxesTest, NonSky2dImageOpens) {
    auto path = GeneratedPath("fdep_time_2d.fits").string();
    auto loader = FileLoader::GetLoader(path);
    loader->OpenFile("0");
    std::string message;
    ASSERT_TRUE(loader->FindCoordinateAxes(message)) << message;
    auto axes = loader->GetAxes();
    EXPECT_EQ(axes.x, 0);
    EXPECT_EQ(axes.y, 1);
    EXPECT_EQ(axes.z, -1);
    Frame frame(0, loader, "0");
    ASSERT_TRUE(frame.IsValid()) << frame.GetErrorMessage();
    CheckPlane(frame, {0, 0});
}

TEST_F(NdAxesTest, CfitsioFallbackReadsAnyRank) {
    // CartaFitsImage is the cfitsio path used for gz/remote/large-header files
    auto path = GeneratedPath("ra_dec_freq_time_beam_stokes.fits").string();
    CartaFitsImage image(path, 0);
    auto shape = image.shape();
    ASSERT_EQ(shape.size(), 6);
    EXPECT_EQ(shape(4), 4);
    EXPECT_EQ(shape(5), 3);
    EXPECT_EQ(image.coordinates().nPixelAxes(), 6);

    casacore::Array<float> data;
    casacore::IPosition start(6, 3, 2, 1, 4, 2, 1), end(6, 3, 2, 1, 4, 2, 1);
    image.doGetSlice(data, casacore::Slicer(start, end, casacore::Slicer::endIsLast));
    EXPECT_FLOAT_EQ(data(casacore::IPosition(6, 0)), EncodedValue({3, 2, 1, 4, 2, 1}));
}

TEST_F(NdAxesTest, FileInfoDescribesAllAxes) {
    auto path = GeneratedPath("ra_dec_freq_stokes_time.fits").string();
    auto loader = FileLoader::GetLoader(path);
    FileExtInfoLoader ext_info_loader(loader);
    std::map<std::string, CARTA::FileInfoExtended> info_map;
    std::string message;
    ASSERT_TRUE(ext_info_loader.FillFitsFileInfoMap(info_map, path, message)) << message;
    ASSERT_EQ(info_map.size(), 1);
    const auto& info = info_map.begin()->second;
    EXPECT_EQ(info.dimensions(), 5);
    EXPECT_EQ(info.width(), 12);
    EXPECT_EQ(info.height(), 10);
    EXPECT_EQ(info.depth(), 6);
    EXPECT_EQ(info.stokes(), 3);
    EXPECT_EQ(info.axes_numbers().depth(), 3);
    EXPECT_EQ(info.axes_numbers().stokes(), 4);

    std::map<std::string, std::string> computed;
    for (int i = 0; i < info.computed_entries_size(); ++i) {
        computed[info.computed_entries(i).name()] = info.computed_entries(i).value();
    }
    EXPECT_EQ(computed["Shape"], "[12, 10, 6, 3, 5] (RA, DEC, FREQ, STOKES, TIME)");
    EXPECT_EQ(computed["Extra axes"], "TIME [5]");
    EXPECT_EQ(computed["Number of channels"], "6");
    EXPECT_EQ(computed["Number of polarizations"], "3");

    // Per-axis descriptions
    ASSERT_EQ(info.axes_size(), 5);
    EXPECT_EQ(info.axes(0).role(), CARTA::AxisRole::AxisX);
    EXPECT_EQ(info.axes(1).role(), CARTA::AxisRole::AxisY);
    EXPECT_EQ(info.axes(2).role(), CARTA::AxisRole::AxisZ);
    EXPECT_EQ(info.axes(3).role(), CARTA::AxisRole::AxisStokes);
    EXPECT_EQ(info.axes(4).role(), CARTA::AxisRole::AxisExtra);
    EXPECT_EQ(info.axes(4).number(), 5);
    EXPECT_EQ(info.axes(4).ctype(), "TIME");
    EXPECT_EQ(info.axes(4).cunit(), "s");
    EXPECT_EQ(info.axes(4).size(), 5);
    ASSERT_EQ(info.axes(4).world_values_size(), 5);
    EXPECT_DOUBLE_EQ(info.axes(4).world_values(3), 3.0);
    EXPECT_EQ(info.axes(0).world_values_size(), 0); // render axes carry no 1D world values
    ASSERT_EQ(info.axes(2).world_values_size(), 6);
    EXPECT_DOUBLE_EQ(info.axes(2).world_values(2), 1.4e9 + 2e6);
}
