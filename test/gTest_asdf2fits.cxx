//
//  gTest_asdf2fits.cxx
//  asdf2fits
//

#include <gtest/gtest.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include <DSTfits/FITSmanager.h>

#include "ASDFmanager.h"
#include "fits2asdf.h"

using namespace DSL;

namespace
{
    std::string fixture(const std::string& name)
    {
        return std::string(ASDF_FIXTURES) + "/" + name + ".asdf";
    }

    class TempFile
    {
    public:
        explicit TempFile(const std::string& extension)
        {
            static int counter = 0;
            const auto* info = testing::UnitTest::GetInstance()->current_test_info();
            std::string base = std::string(info->test_suite_name()) + "_" + info->name();
            for(char& c : base)
                if(c == '/')
                    c = '_';
            path_ = (std::filesystem::temp_directory_path() / (base + "_" + std::to_string(counter++) + extension)).string();
            std::filesystem::remove(path_);
        }
        ~TempFile() { std::filesystem::remove(path_); }
        const std::string& path() const { return path_; }
    private:
        std::string path_;
    };

    std::string readAll(const std::string& path)
    {
        std::ifstream in(path, std::ios::binary);
        std::stringstream text;
        text << in.rdbuf();
        return text.str();
    }

    template<typename T>
    std::shared_ptr<FITSimg<T>> ramp(std::vector<size_t> axes)
    {
        auto image = std::make_shared<FITSimg<T>>(axes);
        std::valarray<T>& pixels = *image->template GetData<T>();
        for(size_t i = 0; i < pixels.size(); i++)
            pixels[i] = static_cast<T>(i % 100 + 1);
        return image;
    }

    template<typename T>
    std::shared_ptr<FITScolumn<T>> column(const std::string& name, dtype type, const std::vector<T>& values, const std::string& unit = "")
    {
        auto c = std::make_shared<FITScolumn<T>>(name, type, unit);
        for(const T& v : values)
            c->push_back(std::any(v));
        return c;
    }
}

// ---------------------------------------------------------------------- Reading existing files

TEST(Metadata, CubeHeaderComesFromTheArrayDescription)
{
    ASDFmanager file(fixture("cube"));
    ASSERT_EQ(file.NumberOfHeader(), 1);

    auto header = file.GetHeader("cube");
    EXPECT_EQ(header->GetInt32ValueForKey("BITPIX"), 8);
    EXPECT_EQ(header->GetInt32ValueForKey("NAXIS"), 3);
    EXPECT_EQ(header->GetUInt64ValueForKey("NAXIS1"), 10u);
    EXPECT_EQ(header->GetUInt64ValueForKey("NAXIS2"), 10u);
    EXPECT_EQ(header->GetUInt64ValueForKey("NAXIS3"), 10u);
    EXPECT_EQ(header->GetValueForKey("EXTNAME"), "cube");
}

TEST(Metadata, HeaderMatchesTheOneOfTheLoadedImage)
{
    ASDFmanager file(fixture("strides"));
    auto header = file.GetHeader("cube_perm");
    auto cube = file.GetImage("cube_perm");

    for(const char* key : {"BITPIX", "NAXIS", "NAXIS1", "NAXIS2", "NAXIS3", "BSCALE", "BZERO"})
        EXPECT_DOUBLE_EQ(header->GetDoubleValueForKey(key), cube->HDU().GetDoubleValueForKey(key)) << key;
    EXPECT_EQ(header->GetValueForKey("EXTNAME"), cube->GetName());
}

TEST(Read, CubePixels)
{
    ASDFmanager file(fixture("cube"));
    auto cube = file.GetPrimary();
    ASSERT_EQ(cube->GetDimension(), 3u);
    EXPECT_EQ(cube->Nelements(), 1000u);
    ASSERT_NE(cube->GetData<uint8_t>(), nullptr);
    for(size_t i : {size_t(0), size_t(1), size_t(255), size_t(999)})
        EXPECT_EQ(cube->UByteValueAtPixel(i), static_cast<uint8_t>(i & 0xFF));
}

TEST(Read, TypedImageThroughDynamicCast)
{
    ASDFmanager file(fixture("cube"));
    EXPECT_NE(std::dynamic_pointer_cast<FITSimg<uint8_t>>(file.GetPrimary()), nullptr);
    EXPECT_EQ(std::dynamic_pointer_cast<FITSimg<float>>(file.GetPrimary()), nullptr);
}

TEST(Read, MultipleHDUPutsPrimaryFirst)
{
    ASDFmanager file(fixture("multiple_hdu"));
    ASSERT_EQ(file.NumberOfHeader(), 4);
    EXPECT_EQ(file.GetPrimaryHeader()->GetInt32ValueForKey("NAXIS"), 0);
    EXPECT_FALSE(file.GetPrimaryHeader()->Exists("EXTNAME"));
    EXPECT_THROW(file.GetPrimary(), FITSexception);

    EXPECT_DOUBLE_EQ(file.GetImage("IMAGE")->at(0), 256.2);
    EXPECT_DOUBLE_EQ(file.GetImage("IMAGE2")->at(0), 1024.44); // big-endian

    auto table = file.GetTable("TABLE");
    ASSERT_EQ(table->ncols(), 2u);
    ASSERT_EQ(table->nrows(), 2u);
    const auto& ids = table->getColumn("ID")->values<std::string>();
    const auto& xs  = std::as_const(*table->getColumn("X")).values<std::string>();
    EXPECT_EQ(ids[1], "                    2");
    EXPECT_EQ(xs[0],  "                  123");
}

TEST(Read, NamesIgnoreCase)
{
    ASDFmanager file(fixture("multiple_hdu"));
    EXPECT_DOUBLE_EQ(file.GetImage("image")->at(0), 256.2);
}

TEST(Read, CompressedBlocks)
{
    ASDFmanager file(fixture("compressed"));
    for(const char* name : {"bzp2", "lz4", "zlib"})
    {
        auto cube = file.GetImage(name);
        ASSERT_EQ(cube->Nelements(), 409600u) << name;
        EXPECT_EQ(cube->UByteValueAtPixel(size_t(0)), 0) << name;
        EXPECT_EQ(cube->UByteValueAtPixel(size_t(409599)), 255) << name;
    }
}

TEST(Read, ByteOrder)
{
    ASDFmanager file(fixture("byteorder"));
    for(const char* name : {"uint16-big", "uint16-little", "uint32-big", "uint32-little", "uint64-big", "uint64-little"})
        for(size_t i = 0; i < 8; i++)
            EXPECT_DOUBLE_EQ(file.GetImage(name)->at(i), static_cast<double>(i)) << name;
}

TEST(Read, StridesAndAxisOrder)
{
    ASDFmanager file(fixture("strides"));

    // numpy shape (3, 4): NAXIS1 = 4, NAXIS2 = 3
    auto fortran = file.GetImage("fortran");
    EXPECT_EQ(fortran->Size(1), 4u);
    EXPECT_EQ(fortran->Size(2), 3u);
    EXPECT_EQ(fortran->LongValueAtPixel(size_t(11)), 11);

    auto perm = file.GetImage("cube_perm");
    EXPECT_EQ(perm->Size(1), 4u);
    EXPECT_EQ(perm->Size(2), 3u);
    EXPECT_EQ(perm->Size(3), 5u);
    EXPECT_EQ(perm->UShortValueAtPixel(size_t(1)), 5);
    EXPECT_EQ(perm->UShortValueAtPixel(size_t(59)), 59);
}

TEST(Read, Offset)
{
    ASDFmanager file(fixture("offset"));
    EXPECT_EQ(file.GetImage("u1")->UByteValueAtPixel(size_t(0)), 16);
    EXPECT_EQ(file.GetImage("i4_2d")->LongValueAtPixel(size_t(0)), 30);
    EXPECT_EQ(file.GetImage("u2_be")->UShortValueAtPixel(size_t(0)), 5);
    EXPECT_EQ(file.GetImage("zlib")->UByteValueAtPixel(size_t(0)), 17);
}

TEST(Read, InlineArray)
{
    ASDFmanager file(fixture("ndarray-inline"));
    EXPECT_DOUBLE_EQ(file.GetImage("explicit")->at(8), 8.);
}

TEST(Read, OtherFixturesDoNotBreakDiscovery)
{
    for(const char* name : {"roman_wcs", "roman_l2_wcs", "value-types", "scalars", "nested", "datatypes"})
        EXPECT_NO_THROW({ ASDFmanager file(fixture(name)); EXPECT_GE(file.NumberOfHeader(), 1) << name; file.GetPrimaryHeader(); }) << name;
}

TEST(Read, GWCSIsNeitherAnHDUNorKeywords)
{
    TempFile path(".asdf");
    {
        std::ofstream out(path.path());
        out << "#ASDF 1.0.0\n#ASDF_STANDARD 1.6.0\n%YAML 1.1\n%TAG ! tag:stsci.edu:asdf/\n--- !core/asdf-1.1.0\n"
            << "sci:\n  telescope: Roman\n  wcs: !<tag:stsci.edu:gwcs/wcs-1.0.0>\n    name: detector\n"
            << "    lookup: !core/ndarray-1.1.0 {data: [1, 2, 3], datatype: int32, shape: [3]}\n"
            << "  data: !core/ndarray-1.1.0 {data: [[1, 2], [3, 4]], datatype: int32, shape: [2, 2]}\n...\n";
    }
    ASDFmanager file(path.path());
    ASSERT_EQ(file.NumberOfHeader(), 1);

    auto header = file.GetPrimaryHeader();
    EXPECT_EQ(header->GetValueForKey("EXTNAME"), "sci.data");
    EXPECT_EQ(header->GetValueForKey("TELESCOPE"), "Roman");
    EXPECT_FALSE(header->Exists("WCS.NAME"));
    EXPECT_EQ(file.GetPrimary()->at(3), 4.);
}

TEST(Read, ScalarsBecomeKeywordsOfThePrimaryHeader)
{
    TempFile path(".asdf");
    {
        std::ofstream out(path.path());
        out << "#ASDF 1.0.0\n#ASDF_STANDARD 1.6.0\n%YAML 1.1\n%TAG ! tag:stsci.edu:asdf/\n--- !core/asdf-1.1.0\n"
            << "object: M31\nflag: true\ncount: 1\nexposure: 1.5\nnothing: null\n"
            << "meta:\n  instrument:\n    name: NISP\n  list: [1, 2]\n...\n";
    }
    ASDFmanager file(path.path());
    ASSERT_EQ(file.NumberOfHeader(), 1);

    auto header = file.GetPrimaryHeader();
    EXPECT_EQ(header->GetValueForKey("OBJECT"), "M31");
    EXPECT_TRUE(header->GetBoolValueForKey("FLAG"));
    EXPECT_EQ(header->GetInt64ValueForKey("COUNT"), 1);
    EXPECT_EQ(header->GetEntry("COUNT")->second.type(), fLongLong);
    EXPECT_DOUBLE_EQ(header->GetDoubleValueForKey("EXPOSURE"), 1.5);
    EXPECT_EQ(header->GetValueForKey("META.INSTRUMENT.NAME"), "NISP");
    EXPECT_FALSE(header->Exists("NOTHING"));
    EXPECT_FALSE(header->Exists("META.LIST"));
}

// ---------------------------------------------------------------------- Errors

TEST(Errors, MissingFile)
{
    EXPECT_THROW(ASDFmanager("/does/not/exist.asdf"), FITSexception);
}

TEST(Errors, NotAnASDFFile)
{
    TempFile path(".asdf");
    { std::ofstream out(path.path()); out << "this is not an ASDF file\n"; }
    EXPECT_THROW(ASDFmanager file(path.path()), FITSexception);
}

TEST(Errors, UnknownHDU)
{
    ASDFmanager file(fixture("multiple_hdu"));
    EXPECT_THROW(file.GetHeader("NOPE"), FITSexception);
    EXPECT_THROW(file.GetImage("NOPE"), FITSexception);
    EXPECT_THROW(file.GetTable("NOPE"), FITSexception);
    EXPECT_THROW(file.GetHeaderAtIndex(0), FITSexception);
    EXPECT_THROW(file.GetHeaderAtIndex(5), FITSexception);
}

TEST(Errors, WrongKindOfHDU)
{
    ASDFmanager file(fixture("multiple_hdu"));
    EXPECT_THROW(file.GetImage("TABLE"), FITSexception);
    EXPECT_THROW(file.GetTable("IMAGE"), FITSexception);
}

TEST(Errors, ClosedFile)
{
    ASDFmanager file(fixture("cube"));
    auto cube = file.GetPrimary();
    file.Close();
    EXPECT_FALSE(file.isOpen());
    EXPECT_THROW(file.GetPrimaryHeader(), FITSexception);
    EXPECT_EQ(cube->UByteValueAtPixel(size_t(3)), 3); // already loaded: still valid
}

TEST(Errors, ExistingFileNeedsReplace)
{
    TempFile path(".asdf");
    auto image = ramp<float>({4, 3});
    fits2asdf(*image, path.path());
    EXPECT_THROW(fits2asdf(*image, path.path()), FITSexception);
    EXPECT_NO_THROW(fits2asdf(*image, path.path(), true));
}

TEST(Errors, VariableLengthColumn)
{
    TempFile path(".asdf");
    FITStable table;
    table.InsertColumn(column<std::vector<int32_t>>("V", tint, {{1, 2}, {1, 2, 3}}));
    EXPECT_THROW(fits2asdf(table, path.path()), FITSexception);
    EXPECT_FALSE(std::filesystem::exists(path.path()));
}

// ---------------------------------------------------------------------- Image round trip

template<typename T> class ImageRoundTrip : public testing::Test {};
using PixelTypes = testing::Types<int8_t, uint8_t, int16_t, uint16_t, int32_t, uint32_t, int64_t, uint64_t, float, double>;
TYPED_TEST_SUITE(ImageRoundTrip, PixelTypes);

TYPED_TEST(ImageRoundTrip, ValuesTypeAndAxes)
{
    using T = TypeParam;
    TempFile path(".asdf");

    auto image = ramp<T>({5, 4, 3});
    (*image->template GetData<T>())[7] = std::numeric_limits<T>::max();
    (*image->template GetData<T>())[8] = std::numeric_limits<T>::lowest();
    fits2asdf(*image, path.path());

    ASDFmanager file(path.path());
    ASSERT_EQ(file.NumberOfHeader(), 1);
    auto back = file.GetPrimary();

    const std::valarray<T>* pixels = back->template GetData<T>();
    ASSERT_NE(pixels, nullptr);
    EXPECT_EQ(back->Size(1), 5u);
    EXPECT_EQ(back->Size(2), 4u);
    EXPECT_EQ(back->Size(3), 3u);
    ASSERT_EQ(pixels->size(), 60u);
    for(size_t i = 0; i < 60; i++)
        EXPECT_EQ((*pixels)[i], (*image->template GetData<T>())[i]) << i;
    EXPECT_EQ(back->GetBitPerPixel(), image->GetBitPerPixel());
    EXPECT_EQ(back->GetEqBitPerPixel(), image->GetEqBitPerPixel());
}

TYPED_TEST(ImageRoundTrip, LoadedImageCanBeSavedAsFITS)
{
    using T = TypeParam;
    TempFile asdf(".asdf");
    TempFile fits(".fits");

    auto image = ramp<T>({6, 2});
    fits2asdf(*image, asdf.path());

    ASDFmanager file(asdf.path());
    file.GetPrimary()->Write(fits.path());

    FITSmanager reader(fits.path());
    auto back = reader.GetPrimary();
    const std::valarray<T>* pixels = back->template GetData<T>();
    ASSERT_NE(pixels, nullptr);
    for(size_t i = 0; i < pixels->size(); i++)
        EXPECT_EQ((*pixels)[i], (*image->template GetData<T>())[i]) << i;
}

TEST(ImageRoundTrip2D, AxisOrderMatchesAstropyConvention)
{
    TempFile path(".asdf");
    auto image = ramp<float>({5, 2});
    fits2asdf(*image, path.path());

    // Only the whitespace of the sequence depends on the libfyaml version
    std::string text = readAll(path.path());
    EXPECT_NE(text.find("tag:astropy.org:astropy/fits/fits-1.2.0"), std::string::npos);
    text.erase(std::remove_if(text.begin(), text.end(), [](unsigned char c){ return std::isspace(c); }), text.end());
    EXPECT_NE(text.find("shape:[2,5]"), std::string::npos);
}

TEST(Header, KeywordTypesAndCommentsSurvive)
{
    TempFile path(".asdf");
    auto image = ramp<float>({2, 2});
    image->HDU().ValueForKey("OBJECT", std::string("M31 'core'"), fChar, "Target name / object");
    image->HDU().ValueForKey("EXPTIME", 1.23456789012345, "Exposure");
    image->HDU().ValueForKey("NCOMBINE", int64_t(-42), "Frames");
    image->HDU().ValueForKey("BIG", uint64_t(18446744073709551615ULL), "Largest");
    image->HDU().ValueForKey("ONE", int32_t(1));
    image->HDU().ValueForKey("ZERO", int32_t(0));
    image->HDU().ValueForKey("DONE", true, "Finished");
    image->HDU().ValueForKey("NOPE", false);
    image->HDU().ValueForKey("HISTORY", std::string(200, 'h'), fChar);
    image->SetName("SCI");
    fits2asdf(*image, path.path());

    ASDFmanager file(path.path());
    auto header = file.GetHeader("SCI");
    EXPECT_EQ(header->GetValueForKey("OBJECT"), "M31 'core'");
    EXPECT_EQ(header->GetEntry("OBJECT")->second.comment(), "Target name / object");
    EXPECT_EQ(header->GetEntry("OBJECT")->second.type(), fChar);
    EXPECT_DOUBLE_EQ(header->GetDoubleValueForKey("EXPTIME"), 1.23456789012345);
    EXPECT_EQ(header->GetInt64ValueForKey("NCOMBINE"), -42);
    EXPECT_EQ(header->GetEntry("NCOMBINE")->second.comment(), "Frames");
    EXPECT_EQ(header->GetUInt64ValueForKey("BIG"), 18446744073709551615ULL);
    EXPECT_EQ(header->GetInt32ValueForKey("ONE"), 1);
    EXPECT_NE(header->GetEntry("ONE")->second.type(), fBool);
    EXPECT_EQ(header->GetInt32ValueForKey("ZERO"), 0);
    EXPECT_NE(header->GetEntry("ZERO")->second.type(), fBool);
    EXPECT_TRUE(header->GetBoolValueForKey("DONE"));
    EXPECT_EQ(header->GetEntry("DONE")->second.type(), fBool);
    EXPECT_FALSE(header->GetBoolValueForKey("NOPE"));
    EXPECT_NE(header->GetValueForKey("HISTORY").find(std::string(72, 'h')), std::string::npos);

    auto cube = file.GetImage("SCI");
    EXPECT_EQ(cube->GetName(), "SCI");
    EXPECT_EQ(cube->HDU().GetValueForKey("OBJECT"), "M31 'core'");
    EXPECT_TRUE(cube->HDU().Exists("HISTORY"));
}

TEST(Masks, RealImageUsesNaN)
{
    TempFile path(".asdf");
    auto image = ramp<float>({4, 2});
    image->MaskPixel(size_t(3));
    fits2asdf(*image, path.path());

    ASDFmanager file(path.path());
    EXPECT_FALSE(file.GetPrimaryHeader()->Exists("BLANK"));
    auto back = file.GetPrimary();
    EXPECT_TRUE(back->isMasked(size_t(3)));
    EXPECT_FALSE(back->isMasked(size_t(2)));
    EXPECT_TRUE(std::isnan(back->FloatValueAtPixel(size_t(3))));
}

TEST(Masks, IntegerImageUsesBlank)
{
    TempFile path(".asdf");
    auto image = ramp<int32_t>({4, 2});
    image->MaskPixel(size_t(3));
    fits2asdf(*image, path.path());

    ASDFmanager file(path.path());
    EXPECT_EQ(file.GetPrimaryHeader()->GetInt64ValueForKey("BLANK"), std::numeric_limits<int32_t>::min());
    auto back = file.GetPrimary();
    EXPECT_TRUE(back->isMasked(size_t(3)));
    EXPECT_FALSE(back->isMasked(size_t(2)));
}

TEST(Masks, UnsignedImageKeepsZeroValid)
{
    TempFile path(".asdf");
    auto image = ramp<uint16_t>({4, 2});
    (*image->GetData<uint16_t>())[0] = 0;
    image->MaskPixel(size_t(3));
    fits2asdf(*image, path.path());

    ASDFmanager file(path.path());
    auto back = file.GetPrimary();
    EXPECT_TRUE(back->isMasked(size_t(3)));
    EXPECT_FALSE(back->isMasked(size_t(0)));
    EXPECT_EQ(back->UShortValueAtPixel(size_t(0)), 0);
}

TEST(Masks, NoBlankWithoutMaskedPixels)
{
    TempFile path(".asdf");
    auto image = ramp<int32_t>({4, 2});
    image->HDU().ValueForKey("BLANK", int32_t(-999));
    fits2asdf(*image, path.path());

    ASDFmanager file(path.path());
    EXPECT_FALSE(file.GetPrimaryHeader()->Exists("BLANK"));
}

// ---------------------------------------------------------------------- Table round trip

TEST(TableRoundTrip, AllColumnKinds)
{
    TempFile path(".asdf");

    FITStable table("CAT");
    table.InsertColumn(column<int8_t>("I8", tsbyte, {-1, 2, 3}));
    table.InsertColumn(column<uint8_t>("U8", tbyte, {1, 2, 255}));
    table.InsertColumn(column<int16_t>("I16", tshort, {-300, 2, 3}, "m"));
    table.InsertColumn(column<uint16_t>("U16", tushort, {1, 2, 65535}));
    table.InsertColumn(column<int32_t>("I32", tint, {-70000, 2, 3}));
    table.InsertColumn(column<uint32_t>("U32", tuint, {1, 2, 4000000000u}));
    table.InsertColumn(column<int64_t>("I64", tlonglong, {-5000000000LL, 2, 3}));
    table.InsertColumn(column<uint64_t>("U64", tulonglong, {1, 2, 18446744073709551615ULL}));
    table.InsertColumn(column<float>("F", tfloat, {1.5f, -2.5f, 3.5f}));
    table.InsertColumn(column<double>("D", tdouble, {1.25, -2.25, 1e300}, "km/s"));
    table.InsertColumn(column<bool>("B", tlogical, {true, false, true}));
    table.InsertColumn(column<std::string>("S", tstring, {"a", "bbbb", ""}));
    table.InsertColumn(column<FITSform::complex>("C", tcplx, {{1.f, 2.f}, {3.f, 4.f}, {5.f, 6.f}}));
    table.InsertColumn(column<FITSform::dblcomplex>("Z", tdbcplx, {{1., 2.}, {3., 4.}, {5., 6.}}));
    table.InsertColumn(column<std::vector<int16_t>>("VI", tshort, {{1, 2}, {3, 4}, {5, 6}}));
    table.InsertColumn(column<std::vector<double>>("VD", tdouble, {{1., 2., 3.}, {4., 5., 6.}, {7., 8., 9.}}));
    table.InsertColumn(column<std::vector<bool>>("VB", tlogical, {{true, false}, {false, true}, {true, true}}));
    table.InsertColumn(column<std::vector<std::string>>("VS", tstring, {{"a", "bc"}, {"def", "g"}, {"", "h"}}));
    table.InsertColumn(column<FITSform::complexVector>("VC", tcplx, {{{1.f, 2.f}, {3.f, 4.f}}, {{5.f, 6.f}, {7.f, 8.f}}, {{9.f, 1.f}, {2.f, 3.f}}}));
    table.HDU().ValueForKey("OBSERVER", std::string("Me"), fChar, "Who");

    fits2asdf(table, path.path());

    ASDFmanager file(path.path());
    ASSERT_EQ(file.NumberOfHeader(), 2);
    EXPECT_EQ(file.GetPrimaryHeader()->GetInt32ValueForKey("NAXIS"), 0);

    auto back = file.GetTable("CAT");
    ASSERT_EQ(back->nrows(), 3u);
    ASSERT_EQ(back->ncols(), table.ncols());
    EXPECT_EQ(back->GetName(), "CAT");
    EXPECT_EQ(back->HDU().GetValueForKey("OBSERVER"), "Me");

    auto same = [&]<typename T>(const std::string& name) {
        const auto& out = std::as_const(*back->getColumn(name)).template values<T>();
        const auto& in  = std::as_const(*table.getColumn(name)).template values<T>();
        ASSERT_EQ(out.size(), in.size()) << name;
        EXPECT_TRUE(out == in) << name;
    };
    same.operator()<int8_t>("I8");
    same.operator()<uint8_t>("U8");
    same.operator()<int16_t>("I16");
    same.operator()<uint16_t>("U16");
    same.operator()<int32_t>("I32");
    same.operator()<uint32_t>("U32");
    same.operator()<int64_t>("I64");
    same.operator()<uint64_t>("U64");
    same.operator()<float>("F");
    same.operator()<double>("D");
    same.operator()<bool>("B");
    same.operator()<std::string>("S");
    same.operator()<FITSform::complex>("C");
    same.operator()<FITSform::dblcomplex>("Z");
    same.operator()<std::vector<int16_t>>("VI");
    same.operator()<std::vector<double>>("VD");
    same.operator()<std::vector<bool>>("VB");
    same.operator()<std::vector<std::string>>("VS");
    same.operator()<FITSform::complexVector>("VC");

    EXPECT_EQ(back->getColumn("D")->getUnit(), "km/s");
    EXPECT_EQ(back->getColumn("I16")->getUnit(), "m");
    EXPECT_EQ(back->getColumn("VD")->getNelem(), 3);

    auto header = file.GetHeader("CAT");
    EXPECT_EQ(header->GetInt32ValueForKey("TFIELDS"), static_cast<int32_t>(table.ncols()));
    EXPECT_EQ(header->GetUInt64ValueForKey("NAXIS2"), 3u);
    EXPECT_EQ(header->GetValueForKey("TTYPE3"), "I16");
    EXPECT_EQ(header->GetValueForKey("TUNIT3"), "m");
}

TEST(TableRoundTrip, EmptyTableKeepsItsColumns)
{
    TempFile path(".asdf");
    FITStable table("EMPTY");
    table.InsertColumn(column<int32_t>("A", tint, {}));
    table.InsertColumn(column<double>("B", tdouble, {}));
    fits2asdf(table, path.path());

    ASDFmanager file(path.path());
    auto back = file.GetTable("EMPTY");
    EXPECT_EQ(back->nrows(), 0u);
    EXPECT_EQ(back->ncols(), 2u);
    EXPECT_EQ(back->getColumn(size_t(2))->getName(), "B");
}

TEST(TableRoundTrip, LoadedTableCanBeSavedAsFITS)
{
    TempFile asdf(".asdf");
    TempFile fits(".fits");

    FITStable table("CAT");
    table.InsertColumn(column<uint16_t>("U16", tushort, {1, 2, 65535}));
    table.InsertColumn(column<double>("D", tdouble, {1.25, -2.25, 1e300}, "km/s"));
    table.InsertColumn(column<std::string>("S", tstring, {"a", "bbbb", ""}));
    table.InsertColumn(column<std::vector<float>>("V", tfloat, {{1.f, 2.f}, {3.f, 4.f}, {5.f, 6.f}}));
    fits2asdf(table, asdf.path());

    ASDFmanager file(asdf.path());
    file.GetTable("CAT")->write(fits.path(), 1);

    FITSmanager reader(fits.path());
    auto back = reader.GetTable("CAT");
    ASSERT_EQ(back->nrows(), 3u);
    EXPECT_EQ(std::as_const(*back->getColumn("U16")).values<uint16_t>()[2], 65535);
    EXPECT_DOUBLE_EQ(std::as_const(*back->getColumn("D")).values<double>()[2], 1e300);
    EXPECT_EQ(back->getColumn("D")->getUnit(), "km/s");
    EXPECT_EQ(std::as_const(*back->getColumn("S")).values<std::string>()[1], "bbbb");
    EXPECT_FLOAT_EQ(std::as_const(*back->getColumn("V")).values<std::vector<float>>()[2][1], 6.f);
}

// ---------------------------------------------------------------------- Whole FITS files

TEST(FitsFile, ConvertedFileMatchesTheFITSOne)
{
    TempFile fitsPath(".fits");
    TempFile asdfPath(".asdf");

    {
        auto primary = ramp<float>({30, 20});
        primary->HDU().ValueForKey("OBJECT", std::string("M31"), fChar, "Target");
        primary->HDU().ValueForKey("CTYPE1", std::string("RA---TAN"), fChar);
        primary->HDU().ValueForKey("CTYPE2", std::string("DEC--TAN"), fChar);
        primary->HDU().ValueForKey("CRPIX1", 15.);
        primary->HDU().ValueForKey("CRPIX2", 10.);
        primary->HDU().ValueForKey("CRVAL1", 10.68);
        primary->HDU().ValueForKey("CRVAL2", 41.27);
        primary->HDU().ValueForKey("CDELT1", -0.001);
        primary->HDU().ValueForKey("CDELT2", 0.001);
        primary->Write(fitsPath.path(), true);

        FITSmanager file(fitsPath.path(), false);
        auto second = ramp<uint16_t>({4, 3, 2});
        second->SetName("U16CUBE");
        file.AppendImage(*second);

        auto table = std::make_shared<FITStable>("CAT");
        table->InsertColumn(column<int32_t>("ID", tint, {1, 2, 3}));
        table->InsertColumn(column<double>("FLUX", tdouble, {0.5, 1.5, 2.5}, "Jy"));
        table->InsertColumn(column<std::string>("NAME", tstring, {"a", "bb", "ccc"}));
        file.InsertTable(table);
        file.Close();
    }

    FITSmanager fits(fitsPath.path());
    fits2asdf(fits, asdfPath.path());
    ASDFmanager asdf(asdfPath.path());

    ASSERT_EQ(asdf.NumberOfHeader(), fits.NumberOfHeader());

    auto fitsPrimary = fits.GetPrimary();
    auto asdfPrimary = asdf.GetPrimary();
    ASSERT_EQ(asdfPrimary->Nelements(), fitsPrimary->Nelements());
    for(size_t i = 0; i < fitsPrimary->Nelements(); i++)
        ASSERT_FLOAT_EQ(asdfPrimary->FloatValueAtPixel(i), fitsPrimary->FloatValueAtPixel(i));
    EXPECT_EQ(asdf.GetPrimaryHeader()->GetValueForKey("OBJECT"), "M31");

    // The WCS is rebuilt from the header
    ASSERT_GT(asdfPrimary->getNumberOfWCS(), 0u);
    const worldCoords a = asdfPrimary->WorldCoordinates(std::vector<size_t>{3, 7});
    const worldCoords b = fitsPrimary->WorldCoordinates(std::vector<size_t>{3, 7});
    ASSERT_EQ(a.size(), b.size());
    for(size_t i = 0; i < a.size(); i++)
        EXPECT_DOUBLE_EQ(a[i], b[i]);

    auto fitsCube = fits.GetImage("U16CUBE");
    auto asdfCube = asdf.GetImage("U16CUBE");
    ASSERT_NE(asdfCube->GetData<uint16_t>(), nullptr);
    for(size_t axis : {1, 2, 3})
        EXPECT_EQ(asdfCube->Size(axis), fitsCube->Size(axis));
    for(size_t i = 0; i < fitsCube->Nelements(); i++)
        ASSERT_EQ(asdfCube->UShortValueAtPixel(i), fitsCube->UShortValueAtPixel(i));

    auto fitsTable = fits.GetTable("CAT");
    auto asdfTable = asdf.GetTable("CAT");
    ASSERT_EQ(asdfTable->nrows(), fitsTable->nrows());
    EXPECT_TRUE(std::as_const(*asdfTable->getColumn("ID")).values<int32_t>() == std::as_const(*fitsTable->getColumn("ID")).values<int32_t>());
    EXPECT_TRUE(std::as_const(*asdfTable->getColumn("FLUX")).values<double>() == std::as_const(*fitsTable->getColumn("FLUX")).values<double>());
    EXPECT_TRUE(std::as_const(*asdfTable->getColumn("NAME")).values<std::string>() == std::as_const(*fitsTable->getColumn("NAME")).values<std::string>());
    EXPECT_EQ(asdfTable->getColumn("FLUX")->getUnit(), "Jy");
}

TEST(FitsFile, ScaledImageIsWrittenAsPhysicalValues)
{
    TempFile fitsPath(".fits");
    TempFile asdfPath(".asdf");
    {
        // BITPIX = 16 with BSCALE and BZERO: the stored integers aren't the pixel values
        fitsfile* raw = nullptr;
        int status = 0;
        long axes[2] = {3, 2};
        short values[6] = {0, 1, 2, 3, 4, 5};
        double scale = 0.5, zero = 10.;
        fits_create_file(&raw, fitsPath.path().c_str(), &status);
        fits_create_img(raw, SHORT_IMG, 2, axes, &status);
        fits_write_key(raw, TDOUBLE, "BSCALE", &scale, "", &status);
        fits_write_key(raw, TDOUBLE, "BZERO", &zero, "", &status);
        fits_set_bscale(raw, 1., 0., &status);
        fits_write_img(raw, TSHORT, 1, 6, values, &status);
        fits_close_file(raw, &status);
        ASSERT_EQ(status, 0);
    }

    FITSmanager fits(fitsPath.path());
    auto original = fits.GetPrimary();
    EXPECT_DOUBLE_EQ(original->at(4), 12.);

    fits2asdf(fits, asdfPath.path());
    ASDFmanager asdf(asdfPath.path());
    auto back = asdf.GetPrimary();
    for(size_t i = 0; i < 6; i++)
        EXPECT_DOUBLE_EQ(back->at(i), original->at(i));
    EXPECT_DOUBLE_EQ(back->HDU().GetDoubleValueForKey("BSCALE"), 1.);
}

// libasdf looped forever when the "\n..." that ends the tree lay at the end of its read buffer: try
// every position of the marker within a buffer (see cmake/patches/). A failure shows up as a timeout.
TEST(ASDFparser, DocumentEndMarkerAtAnyBufferPosition)
{
    const std::string head = "#ASDF 1.0.0\n#ASDF_STANDARD 1.6.0\n%YAML 1.1\n%TAG ! tag:stsci.edu:asdf/\n"
                             "--- !core/asdf-1.1.0\nOBJECT: ";
    TempFile file(".asdf");
    for(size_t pad = 1; pad <= BUFSIZ; pad++)
    {
        {
            std::ofstream out(file.path(), std::ios::binary);
            out << head << std::string(pad, 'x') << "\n...\n" << std::string(64, '\n');
        }
        ASDFmanager asdf(file.path());
        ASSERT_TRUE(asdf.GetPrimaryHeader()->Exists("OBJECT")) << "pad=" << pad;
    }
}

TEST(ASDFparser, DocumentEndMarkerFollowedByOtherBytes)
{
    TempFile file(".asdf");
    {
        std::ofstream out(file.path(), std::ios::binary);
        out << "#ASDF 1.0.0\n#ASDF_STANDARD 1.6.0\n%YAML 1.1\n%TAG ! tag:stsci.edu:asdf/\n"
               "--- !core/asdf-1.1.0\nOBJECT: a\n... not a marker\n...\n";
    }
    // The content isn't valid YAML, so an error is fine: only returning matters.
    try { ASDFmanager asdf(file.path()); asdf.NumberOfHeader(); } catch(const FITSexception&) {}
}
