//
//  gTest_asdf2fits.cxx
//  asdf2fits
//

#include <gtest/gtest.h>

#include <algorithm>
#include <bit>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <filesystem>
#include <fstream>
#include <functional>
#include <limits>
#include <locale>
#include <memory>
#include <numbers>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#include <cstring>

#include <asdf.h>

#include <DSTfits/FITSmanager.h>

#include "ASDFcube.h"
#include "ASDFmanager.h"
#include "fits2asdf.h"

using namespace DSL;

namespace
{
    std::string fixture(const std::string& name)
    {
        return std::string(ASDF_FIXTURES) + "/" + name + ".asdf";
    }

    // Files written by the tests go to build/Testing/Temporary (ASDF2FITS_TEST_DIR, set by CMakeLists.txt), never
    // to the temporary directory of the system. They are removed at the end of each test, unless
    // ASDF2FITS_KEEP_TEST_FILES is set:
    //   unset, empty or 0  removed;
    //   1                  kept in build/Testing/Temporary;
    //   anything else      the directory to keep them in, created if needed (relative to the working directory).
    // Each file kept is printed when its test ends.
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
            path_ = (directory() / (base + "_" + std::to_string(counter++) + extension)).string();
            std::filesystem::remove(path_);
        }
        ~TempFile()
        {
            if(!keep())
                std::filesystem::remove(path_);
            else if(std::filesystem::exists(path_))
                std::cout << "[   KEPT   ] " << path_ << std::endl;
        }
        const std::string& path() const { return path_; }

    private:
        static const char* setting()
        {
            const char* value = std::getenv("ASDF2FITS_KEEP_TEST_FILES");
            return (value && *value && std::string(value) != "0") ? value : nullptr;
        }
        static bool keep() { return setting() != nullptr; }
        static std::filesystem::path directory()
        {
            const char* value = setting();
            const std::filesystem::path chosen = (!value || std::string(value) == "1") ? std::filesystem::path(ASDF2FITS_TEST_DIR)
                                                                                       : std::filesystem::path(value);
            std::error_code error;
            if(!std::filesystem::is_directory(chosen, error))
                std::filesystem::create_directories(chosen, error);
            if(!std::filesystem::is_directory(chosen, error))
                throw std::runtime_error("can't create the directory of the test files: " + chosen.string());
            return chosen;
        }

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

    // ------------------------------------------------------------------ Files written with libasdf alone

    struct TestBlock
    {
        std::vector<unsigned char> bytes;
        std::string compression;
    };

    template<typename T>
    TestBlock blockOf(const std::vector<T>& values, const std::string& compression = "")
    {
        TestBlock block;
        block.bytes.resize(values.size() * sizeof(T));
        std::memcpy(block.bytes.data(), values.data(), block.bytes.size());
        block.compression = compression;
        return block;
    }

    // The tree is the YAML below the "--- !core/asdf-1.1.0" line; blocks are little-endian
    void writeTestFile(const std::string& path, const std::string& tree, const std::vector<TestBlock>& blocks)
    {
        const std::string yaml = "#ASDF 1.0.0\n#ASDF_STANDARD 1.6.0\n%YAML 1.1\n%TAG ! tag:stsci.edu:asdf/\n--- !core/asdf-1.1.0\n" + tree + "...\n";
        asdf_file_t* file = asdf_open_mem(yaml.data(), yaml.size());
        ASSERT_NE(file, nullptr);
        for(const TestBlock& block : blocks)
        {
            asdf_block_t* b = asdf_block_create(block.bytes.data(), block.bytes.size());
            ASSERT_NE(b, nullptr);
            if(!block.compression.empty())
                ASSERT_EQ(asdf_block_compression_set(b, block.compression.c_str()), 0);
            ASSERT_NE(asdf_block_append(file, b), nullptr);
        }
        ASSERT_EQ(asdf_write_to_file(file, path.c_str()), 0);
        asdf_close(file);
    }

    // A file laid out as a Roman data model: an image with err, dq and a calibration array of another shape
    constexpr size_t kNX = 64;
    constexpr size_t kNY = 48;
    float    dataAt (size_t x, size_t y) { return static_cast<float>(y * 1000 + x); }
    float    errAt  (size_t x, size_t y) { return static_cast<float>(x + y) + 0.25f; }
    uint32_t dqAt   (size_t x, size_t y) { return (x % 7 == 3 ? 1u : 0u) | (y == 5 ? 4u : 0u); }
    uint16_t amp33At(size_t x, size_t y) { return static_cast<uint16_t>(y * 10 + x); }

    // A GWCS of one frame and no transform: it can't be used, so changes of shape remove it
    const char* const kStubWcs =
        "    wcs: !<tag:stsci.edu:gwcs/wcs-1.4.0>\n"
        "      name: detector\n"
        "      steps:\n"
        "      - !<tag:stsci.edu:gwcs/step-1.3.0>\n"
        "        frame: !<tag:stsci.edu:gwcs/frame2d-1.2.0>\n"
        "          name: detector\n"
        "          unit: [!unit/unit-1.0.0 pixel, !unit/unit-1.0.0 pixel]\n";

    // A linear GWCS of the 64x48 image: u = (x - 31.5) / 1000, v = (y - 23.5) / 500, in the bounding box
    const char* const kLinearWcs =
        "    wcs: !<tag:stsci.edu:gwcs/wcs-1.4.0>\n"
        "      name: linear\n"
        "      pixel_shape: null\n"
        "      steps:\n"
        "      - !<tag:stsci.edu:gwcs/step-1.3.0>\n"
        "        frame: !<tag:stsci.edu:gwcs/frame2d-1.2.0>\n"
        "          axes_names: [x, y]\n"
        "          name: detector\n"
        "          unit: [!unit/unit-1.0.0 pixel, !unit/unit-1.0.0 pixel]\n"
        "        transform: !transform/compose-1.3.0\n"
        "          bounding_box: !transform/property/bounding_box-1.1.0\n"
        "            ignore: []\n"
        "            intervals:\n"
        "              x0: [-0.5, 63.5]\n"
        "              x1: [-0.5, 47.5]\n"
        "            order: F\n"
        "          forward:\n"
        "          - !transform/concatenate-1.3.0\n"
        "            forward:\n"
        "            - !transform/shift-1.3.0 {offset: -31.5}\n"
        "            - !transform/shift-1.3.0 {offset: -23.5}\n"
        "          - !transform/concatenate-1.3.0\n"
        "            forward:\n"
        "            - !transform/scale-1.3.0 {factor: 0.001}\n"
        "            - !transform/scale-1.3.0 {factor: 0.002}\n"
        "          inputs: [x0, x1]\n"
        "          outputs: [y0, y1]\n"
        "      - !<tag:stsci.edu:gwcs/step-1.3.0>\n"
        "        frame: !<tag:stsci.edu:gwcs/frame2d-1.2.0>\n"
        "          axes_names: [u, v]\n"
        "          name: world\n"
        "          unit: [!unit/unit-1.0.0 deg, !unit/unit-1.0.0 deg]\n"
        "        transform: null\n";

    // u = 1 + x/2 + y/1000 + x^2/10^4, v = -2 + x/500 + 0.6 y + x y/5000: polynomials without an inverse
    const char* const kPolynomialWcs =
        "    wcs: !<tag:stsci.edu:gwcs/wcs-1.4.0>\n"
        "      name: polynomial\n"
        "      steps:\n"
        "      - !<tag:stsci.edu:gwcs/step-1.3.0>\n"
        "        frame: !<tag:stsci.edu:gwcs/frame2d-1.2.0> {name: detector}\n"
        "        transform: !transform/compose-1.3.0\n"
        "          forward:\n"
        "          - !transform/remap_axes-1.4.0 {mapping: [0, 1, 0, 1]}\n"
        "          - !transform/concatenate-1.3.0\n"
        "            forward:\n"
        "            - !transform/polynomial-1.2.0\n"
        "              coefficients: [[1.0, 0.001, 0.0], [0.5, 0.0, 0.0], [0.0001, 0.0, 0.0]]\n"
        "            - !transform/polynomial-1.2.0\n"
        "              coefficients: [[-2.0, 0.6, 0.0], [0.002, 0.0002, 0.0], [0.0, 0.0, 0.0]]\n"
        "      - !<tag:stsci.edu:gwcs/step-1.3.0>\n"
        "        frame: !<tag:stsci.edu:gwcs/frame2d-1.2.0> {name: world}\n"
        "        transform: null\n";

    void writeRoman(const std::string& path, const std::string& wcs = kStubWcs)
    {
        std::vector<float> data, err;
        std::vector<uint32_t> dq;
        std::vector<uint16_t> amp33;
        for(size_t y = 0; y < kNY; y++)
        {
            for(size_t x = 0; x < kNX; x++)
            {
                data.push_back(dataAt(x, y));
                err.push_back(errAt(x, y));
                dq.push_back(dqAt(x, y));
            }
            for(size_t x = 0; x < 4; x++)
                amp33.push_back(amp33At(x, y));
        }

        const std::string tree =
            "history:\n"
            "  entries:\n"
            "  - !core/history_entry-1.0.0\n"
            "    description: made by the test\n"
            "    time: '2020-01-01T00:00:00'\n"
            "    software: !core/software-1.0.0 {name: test, version: '1.0'}\n"
            "roman: !<asdf://stsci.edu/datamodels/roman/tags/wfi_image-1.4.0>\n"
            "  meta:\n"
            "    exposure:\n"
            "      start_time: !time/time-1.4.0 2027-06-01T00:00:00.000\n"
            "      exposure_time: 133.76\n"
            "      nresultants: 10\n"
            "    instrument: {detector: WFI01, name: WFI, optical_element: F158}\n"
            "    cal_logs: ['first line', \"second line\\nwith a break\"]\n"
            + wcs +
            "  data: !core/ndarray-1.1.0 {source: 0, datatype: float32, byteorder: little, shape: [48, 64]}\n"
            "  err: !core/ndarray-1.1.0 {source: 1, datatype: float32, byteorder: little, shape: [48, 64]}\n"
            "  dq: !core/ndarray-1.1.0 {source: 2, datatype: uint32, byteorder: little, shape: [48, 64]}\n"
            "  amp33: !core/ndarray-1.1.0 {source: 3, datatype: uint16, byteorder: little, shape: [48, 4]}\n";
        writeTestFile(path, tree, {blockOf(data), blockOf(err, "lz4"), blockOf(dq, "zlib"), blockOf(amp33, "bzp2")});
    }

    // The compression and the bytes, as stored, of the block of an ndarray
    std::pair<std::string, std::string> rawBlock(const std::string& path, const std::string& array)
    {
        asdf_file_t* file = asdf_open_file(path.c_str(), "r");
        if(!file)
            return {};
        asdf_ndarray_t* ndarray = nullptr;
        std::pair<std::string, std::string> result;
        if(asdf_get_ndarray(file, array.c_str(), &ndarray) == ASDF_VALUE_OK)
        {
            asdf_block_t* block = asdf_block_open(file, ndarray->source);
            size_t size = 0;
            const void* raw = asdf_block_data_raw(block, &size);
            result.first = asdf_block_compression(block);
            result.second.assign(static_cast<const char*>(raw), size);
            asdf_block_close(block);
            asdf_ndarray_destroy(ndarray);
        }
        asdf_close(file);
        return result;
    }

    // Two trees hold the same: tags, keys, values, and arrays of the same shape, type and pixels
    void expectSameTree(const ASDFnode& a, const ASDFnode& b, const std::string& path = "")
    {
        if(a.IsArray() || b.IsArray())
        {
            ASSERT_TRUE(a.IsArray() && b.IsArray()) << path;
            EXPECT_EQ(a.Array()->Shape(), b.Array()->Shape()) << path;
            EXPECT_EQ(a.Array()->DataType(), b.Array()->DataType()) << path;
            if(a.Array()->IsImage())
            {
                auto pa = a.Array()->Read({});
                auto pb = b.Array()->Read({});
                ASSERT_EQ(pa->Nelements(), pb->Nelements()) << path;
                for(size_t i = 0; i < pa->Nelements(); i++)
                    ASSERT_TRUE(pa->at(i) == pb->at(i) || (std::isnan(pa->at(i)) && std::isnan(pb->at(i)))) << path << " [" << i << "]";
            }
            return;
        }
        EXPECT_EQ(a.Tag(), b.Tag()) << path;
        ASSERT_EQ(a.GetKind(), b.GetKind()) << path;
        if(a.IsMapping() || a.IsSequence())
        {
            ASSERT_EQ(a.Keys(), b.Keys()) << path;
            ASSERT_EQ(a.Size(), b.Size()) << path;
            for(size_t i = 0; i < a.Size(); i++)
                expectSameTree(a.Items()[i], b.Items()[i], path + "/" + (a.IsMapping() ? a.Keys()[i] : std::to_string(i)));
            return;
        }
        EXPECT_TRUE(a == b) << path << ": " << a.Text() << " vs " << b.Text();
    }

    std::shared_ptr<ASDFcube> asASDF(const std::shared_ptr<FITScube>& cube)
    {
        return std::dynamic_pointer_cast<ASDFcube>(cube);
    }

    // ------------------------------------------------------------------ Multi-HDU file built with the libasdf API
    //
    // Every pixel type has its own HDU, of its own shape, byte order and compression, so that one file covers
    // what a reader has to get right. Nothing below goes through asdf2fits: the arrays are written with
    // asdf_set_ndarray() and read back with asdf_ndarray_read_all().

    template<typename F>
    void forEachPixelType(F&& f)
    {
        f.template operator()<int8_t>();
        f.template operator()<uint8_t>();
        f.template operator()<int16_t>();
        f.template operator()<uint16_t>();
        f.template operator()<int32_t>();
        f.template operator()<uint32_t>();
        f.template operator()<int64_t>();
        f.template operator()<uint64_t>();
        f.template operator()<float>();
        f.template operator()<double>();
    }

    template<typename T>
    constexpr asdf_scalar_datatype_t asdfTypeOf()
    {
        if constexpr(std::is_same_v<T, int8_t>)   return ASDF_DATATYPE_INT8;
        else if constexpr(std::is_same_v<T, uint8_t>)  return ASDF_DATATYPE_UINT8;
        else if constexpr(std::is_same_v<T, int16_t>)  return ASDF_DATATYPE_INT16;
        else if constexpr(std::is_same_v<T, uint16_t>) return ASDF_DATATYPE_UINT16;
        else if constexpr(std::is_same_v<T, int32_t>)  return ASDF_DATATYPE_INT32;
        else if constexpr(std::is_same_v<T, uint32_t>) return ASDF_DATATYPE_UINT32;
        else if constexpr(std::is_same_v<T, int64_t>)  return ASDF_DATATYPE_INT64;
        else if constexpr(std::is_same_v<T, uint64_t>) return ASDF_DATATYPE_UINT64;
        else if constexpr(std::is_same_v<T, float>)    return ASDF_DATATYPE_FLOAT32;
        else if constexpr(std::is_same_v<T, double>)   return ASDF_DATATYPE_FLOAT64;
        else static_assert(sizeof(T) == 0, "not a pixel type");
    }

    template<typename T>
    std::string asdfNameOf()
    {
        if constexpr(std::is_same_v<T, float>)  return "float32";
        else if constexpr(std::is_same_v<T, double>) return "float64";
        else
            return std::string(std::is_signed_v<T> ? "int" : "uint") + std::to_string(8 * sizeof(T));
    }

    // One HDU of the file: the ndarray is at <group>/data, next to the scalars <group>/label and <group>/rank
    struct ArraySpec
    {
        std::string group;
        std::vector<uint64_t> shape;    // numpy order: the last one is NAXIS1
        asdf_byteorder_t order;
        std::string compression;        // "" for none
        uint64_t salt;                  // makes the pixels of each HDU different
    };

    template<typename T>
    ArraySpec specOf()
    {
        if constexpr(std::is_same_v<T, int8_t>)   return {"i8",  {6},          ASDF_BYTEORDER_LITTLE, "",     1};
        else if constexpr(std::is_same_v<T, uint8_t>)  return {"u8",  {5, 7},       ASDF_BYTEORDER_LITTLE, "bzp2", 2};
        else if constexpr(std::is_same_v<T, int16_t>)  return {"i16", {3, 4, 5},    ASDF_BYTEORDER_BIG,    "",     3};
        else if constexpr(std::is_same_v<T, uint16_t>) return {"u16", {7, 6},       ASDF_BYTEORDER_LITTLE, "lz4",  4};
        else if constexpr(std::is_same_v<T, int32_t>)  return {"i32", {2, 3, 4, 5}, ASDF_BYTEORDER_LITTLE, "zlib", 5};
        else if constexpr(std::is_same_v<T, uint32_t>) return {"u32", {9, 4},       ASDF_BYTEORDER_BIG,    "",     6};
        else if constexpr(std::is_same_v<T, int64_t>)  return {"i64", {3, 2, 4},    ASDF_BYTEORDER_LITTLE, "",     7};
        else if constexpr(std::is_same_v<T, uint64_t>) return {"u64", {11},         ASDF_BYTEORDER_BIG,    "zlib", 8};
        else if constexpr(std::is_same_v<T, float>)    return {"f32", {8, 5},       ASDF_BYTEORDER_BIG,    "lz4",  9};
        else if constexpr(std::is_same_v<T, double>)   return {"f64", {4, 3, 2},    ASDF_BYTEORDER_LITTLE, "bzp2", 10};
        else static_assert(sizeof(T) == 0, "not a pixel type");
    }

    std::vector<size_t> axesOf(const std::vector<uint64_t>& shape)
    {
        return std::vector<size_t>(shape.rbegin(), shape.rend());
    }

    // The extremes of the type first, then values of both signs that differ from one pixel to the next
    template<typename T>
    T pixelValue(size_t i, uint64_t salt)
    {
        if(i == 0)
            return std::numeric_limits<T>::lowest();
        if(i == 1)
            return std::numeric_limits<T>::max();
        const uint64_t magnitude = 3 * i + salt;
        if constexpr(std::is_floating_point_v<T>)
            return static_cast<T>((i % 2 ? -0.375 : 0.625) * static_cast<double>(magnitude));
        else if constexpr(std::is_signed_v<T>)
            return static_cast<T>(i % 2 ? -static_cast<int64_t>(magnitude) : static_cast<int64_t>(magnitude));
        else
            return static_cast<T>(magnitude);
    }

    // Pixels of an HDU, in file order (the last axis of the shape varies fastest), which is also FITS order
    template<typename T>
    std::vector<T> pixelsOf(const ArraySpec& spec, uint64_t salt)
    {
        size_t n = 1;
        for(uint64_t length : spec.shape)
            n *= length;
        std::vector<T> values(n);
        for(size_t i = 0; i < n; i++)
            values[i] = pixelValue<T>(i, salt);
        return values;
    }

    template<typename T>
    std::vector<T> pixelsOf(const ArraySpec& spec) { return pixelsOf<T>(spec, spec.salt); }

    // libasdf 0.3.0 can't add a key to a mapping that is already in the tree: asdf_set_*(file, "a/b") fails
    // with ASDF_VALUE_ERR_OOM when "a" exists (asdf_node_materialize_path appends "a" to its parent again).
    // So each mapping is built whole, as asdf2fits does, and put in the tree once.
    // LibasdfKnownBugs.SetUnderAnExistingMapping fails once a libasdf without the bug is built.
    template<typename T>
    void setNdarray(asdf_file_t* file, asdf_mapping_t* mapping, const std::string& key, const std::vector<uint64_t>& shape,
                    const std::vector<T>& values, asdf_byteorder_t order, const std::string& compression)
    {
        asdf_ndarray_t ndarray{};
        ndarray.datatype.type = asdfTypeOf<T>();
        ndarray.byteorder = order;
        ndarray.ndim = static_cast<uint32_t>(shape.size());
        ndarray.shape = shape.data();

        // The bytes are given as stored: in the byte order announced
        std::vector<unsigned char> bytes(values.size() * sizeof(T));
        std::memcpy(bytes.data(), values.data(), bytes.size());
        if((order == ASDF_BYTEORDER_BIG) != (std::endian::native == std::endian::big))
            for(size_t i = 0; i < bytes.size(); i += sizeof(T))
                std::reverse(bytes.begin() + static_cast<std::ptrdiff_t>(i), bytes.begin() + static_cast<std::ptrdiff_t>(i + sizeof(T)));

        ASSERT_EQ(asdf_ndarray_data_copy(&ndarray, bytes.data()), ASDF_NDARRAY_OK) << key;
        if(!compression.empty() && asdf_ndarray_compression_set(&ndarray, compression.c_str()) != 0)
        {
            asdf_ndarray_data_dealloc(&ndarray);
            FAIL() << key << ": compression " << compression << " refused";
        }
        // The file takes the data: it is freed when the file is closed
        asdf_value_t* value = asdf_value_of_ndarray(file, &ndarray);
        if(!value)
        {
            asdf_ndarray_data_dealloc(&ndarray);
            FAIL() << key << ": " << (asdf_error(file) ? asdf_error(file) : "asdf_value_of_ndarray failed");
        }
        ASSERT_EQ(asdf_mapping_set(mapping, key.c_str(), value), ASDF_VALUE_OK) << key;
    }

    // Put a mapping built with asdf_mapping_create() at a key of the root; the file takes it
    void setGroup(asdf_file_t* file, const std::string& key, asdf_mapping_t* mapping)
    {
        ASSERT_NE(mapping, nullptr) << key;
        ASSERT_EQ(asdf_set_mapping(file, key.c_str(), mapping), ASDF_VALUE_OK) << key;
    }

    // A mapping whose image has companions, as in the data models of Roman or JWST
    constexpr size_t kSciNX = 8;
    constexpr size_t kSciNY = 6;
    float    sciDataAt(size_t x, size_t y) { return static_cast<float>(y * 100 + x) + 0.5f; }
    float    sciErrAt (size_t x, size_t y) { return 0.25f * static_cast<float>(x + y) + 1.f; }
    uint16_t sciDqAt  (size_t x, size_t y) { return static_cast<uint16_t>((x == 2 ? 1 : 0) | (y == 4 ? 2 : 0)); }

    constexpr int kMultiHDUs = 10 + 3;      // one per pixel type, and sci/data, sci/err, sci/dq

    void writeMultiHDU(const std::string& path)
    {
        asdf_file_t* file = asdf_open_mem(nullptr, 0);
        ASSERT_NE(file, nullptr);

        // Scalars at the root: keywords of the primary HDU
        EXPECT_EQ(asdf_set_string0(file, "object", "NGC 1300"), ASDF_VALUE_OK);
        EXPECT_EQ(asdf_set_double(file, "exptime", 87.5), ASDF_VALUE_OK);
        EXPECT_EQ(asdf_set_int64(file, "ncombine", -12), ASDF_VALUE_OK);
        EXPECT_EQ(asdf_set_bool(file, "calibrated", true), ASDF_VALUE_OK);

        forEachPixelType([&]<typename T>()
        {
            const ArraySpec spec = specOf<T>();
            asdf_mapping_t* group = asdf_mapping_create(file);
            ASSERT_NE(group, nullptr) << spec.group;
            EXPECT_EQ(asdf_mapping_set_string0(group, "label", (spec.group + " image").c_str()), ASDF_VALUE_OK);
            EXPECT_EQ(asdf_mapping_set_int64(group, "rank", static_cast<int64_t>(spec.shape.size())), ASDF_VALUE_OK);
            setNdarray<T>(file, group, "data", spec.shape, pixelsOf<T>(spec), spec.order, spec.compression);
            setGroup(file, spec.group, group);
        });

        std::vector<float> data, err;
        std::vector<uint16_t> dq;
        for(size_t y = 0; y < kSciNY; y++)
            for(size_t x = 0; x < kSciNX; x++)
            {
                data.push_back(sciDataAt(x, y));
                err.push_back(sciErrAt(x, y));
                dq.push_back(sciDqAt(x, y));
            }
        asdf_mapping_t* sci = asdf_mapping_create(file);
        ASSERT_NE(sci, nullptr);
        EXPECT_EQ(asdf_mapping_set_string0(sci, "telescope", "Euclid"), ASDF_VALUE_OK);
        setNdarray<float>(file, sci, "data", {kSciNY, kSciNX}, data, ASDF_BYTEORDER_LITTLE, "");
        setNdarray<float>(file, sci, "err", {kSciNY, kSciNX}, err, ASDF_BYTEORDER_BIG, "lz4");
        setNdarray<uint16_t>(file, sci, "dq", {kSciNY, kSciNX}, dq, ASDF_BYTEORDER_LITTLE, "zlib");
        setGroup(file, "sci", sci);

        EXPECT_EQ(asdf_write_to_file(file, path.c_str()), 0);
        asdf_close(file);
    }

    // An ndarray as libasdf alone reads it: pixels in the native byte order, of the type stored
    struct LibasdfArray
    {
        asdf_scalar_datatype_t type = ASDF_DATATYPE_UNKNOWN;
        std::vector<uint64_t> shape;
        std::string compression;
        std::vector<unsigned char> bytes;

        template<typename T>
        std::vector<T> values() const
        {
            std::vector<T> out(bytes.size() / sizeof(T));
            std::memcpy(out.data(), bytes.data(), out.size() * sizeof(T));
            return out;
        }
    };

    void readNdarray(const std::string& path, const std::string& array, LibasdfArray& out)
    {
        asdf_file_t* file = asdf_open_file(path.c_str(), "r");
        ASSERT_NE(file, nullptr) << path;
        asdf_ndarray_t* ndarray = nullptr;
        if(asdf_get_ndarray(file, array.c_str(), &ndarray) != ASDF_VALUE_OK || !ndarray)
        {
            asdf_close(file);
            FAIL() << "no ndarray at " << array << " in " << path;
        }
        out.type = ndarray->datatype.type;
        out.shape.assign(ndarray->shape, ndarray->shape + ndarray->ndim);
        if(asdf_block_t* block = asdf_block_open(file, ndarray->source))
        {
            out.compression = asdf_block_compression(block);
            asdf_block_close(block);
        }
        void* data = nullptr;
        const asdf_ndarray_err_t status = asdf_ndarray_read_all(ndarray, ASDF_DATATYPE_SOURCE, &data);
        if(status == ASDF_NDARRAY_OK && data)
        {
            const auto* first = static_cast<const unsigned char*>(data);
            out.bytes.assign(first, first + asdf_ndarray_nbytes(ndarray));
        }
        asdf_free(data);
        asdf_ndarray_destroy(ndarray);
        asdf_close(file);
        ASSERT_EQ(status, ASDF_NDARRAY_OK) << array;
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

// ---------------------------------------------------------------------- ASDF trees

TEST(ASDFnode, PathsKindsAndValues)
{
    const std::string key = "x/y~z";
    const std::string path = ASDFnode::ChildPath("a", key);
    EXPECT_EQ(path, "a/x~1y~0z");
    EXPECT_EQ(ASDFnode::SplitPath(path), (std::vector<std::string>{"a", key}));
    EXPECT_EQ(ASDFnode::ParentPath(path), "a");

    ASDFnode tree;
    tree.SetPath(path, 3);
    ASSERT_TRUE(tree.IsMapping());
    EXPECT_EQ(tree["a"][key].Get<int>(), 3);
    EXPECT_EQ(tree.At(path).Get<double>(), 3.);
    EXPECT_THROW(tree.At("a/missing"), FITSexception);
    EXPECT_THROW(tree.At(path).Get<std::string>(), FITSexception);
    EXPECT_THROW(ASDFnode(300).Get<int8_t>(), FITSexception);
    EXPECT_TRUE(ASDFnode(int64_t(5)) == ASDFnode(uint64_t(5)));
    EXPECT_FALSE(ASDFnode(5) == ASDFnode(5.));

    // A time replaced by another scalar stays a time
    tree.SetPath("when", ASDFnode("2020-01-01"));
    tree.At("when").SetTag("tag:stsci.edu:asdf/time/time-1.4.0");
    tree.SetPath("when", ASDFnode("2021-01-01"));
    EXPECT_EQ(tree.At("when").Tag(), "tag:stsci.edu:asdf/time/time-1.4.0");
    EXPECT_EQ(tree.At("when").Get<std::string>(), "2021-01-01");

    ASDFnode& list = tree["list"];
    list.Append(1);
    list.Append("two");
    EXPECT_EQ(tree.At("list/1").Get<std::string>(), "two");
    EXPECT_EQ(tree.Keys(), (std::vector<std::string>{"a", "when", "list"}));

    EXPECT_TRUE(tree.ErasePath("list/0"));
    EXPECT_EQ(tree.At("list/0").Get<std::string>(), "two");
    EXPECT_TRUE(tree.ErasePath(path));
    EXPECT_EQ(tree.FindPath(path), nullptr);

    // Reals are written so that they read back as reals
    EXPECT_EQ(ASDFnode(2.).Text(), "2.0");
    EXPECT_EQ(ASDFnode(1e20).Text(), "1.0e+20");
    EXPECT_EQ(ASDFnode(std::nan("")).Text(), ".nan");
}

TEST(ASDFtree, ReadKeepsTagsAndTypesAndReadsNoArray)
{
    TempFile path(".asdf");
    writeRoman(path.path());
    ASDFmanager file(path.path());
    ASDFnode tree = file.GetTree();

    EXPECT_EQ(tree.At("roman").Tag(), "asdf://stsci.edu/datamodels/roman/tags/wfi_image-1.4.0");
    const ASDFnode& start = tree.At("roman/meta/exposure/start_time");
    EXPECT_EQ(start.Get<std::string>(), "2027-06-01T00:00:00.000");
    EXPECT_EQ(start.Tag(), "tag:stsci.edu:asdf/time/time-1.4.0");
    EXPECT_EQ(tree.At("roman/meta/exposure/nresultants").Get<int>(), 10);
    EXPECT_DOUBLE_EQ(tree.At("roman/meta/exposure/exposure_time").Get<double>(), 133.76);
    EXPECT_EQ(tree.At("roman/meta/cal_logs/1").Get<std::string>(), "second line\nwith a break");
    EXPECT_TRUE(tree.At("roman/meta/wcs").HasTag("gwcs/wcs-"));
    EXPECT_EQ(tree.At("roman/meta/wcs/steps/0/frame/unit/1").Tag(), "tag:stsci.edu:asdf/unit/unit-1.0.0");

    auto dq = tree.At("roman/dq").Array();
    ASSERT_TRUE(dq);
    EXPECT_FALSE(dq->IsLoaded());
    EXPECT_TRUE(dq->IsInFile());
    EXPECT_TRUE(dq->IsImage());
    EXPECT_EQ(dq->Shape(), (std::vector<uint64_t>{kNY, kNX}));
    EXPECT_EQ(dq->Axes(), (std::vector<size_t>{kNX, kNY}));
    EXPECT_EQ(dq->DataType(), "uint32");
    EXPECT_EQ(dq->Bytes(), kNX * kNY * 4);
    EXPECT_EQ(dq->Compression(), "zlib");
    EXPECT_EQ(tree.At("roman/data").Array()->Compression(), "");
    EXPECT_DOUBLE_EQ(dq->ValueAt({3, 5}), dqAt(3, 5));
    EXPECT_FALSE(dq->IsLoaded());

    auto header = file.GetHeader("roman.data");
    EXPECT_EQ(header->GetInt64ValueForKey("META.EXPOSURE.NRESULTANTS"), 10);
    EXPECT_EQ(header->GetValueForKey("META.INSTRUMENT.NAME"), "WFI");
}

// ---------------------------------------------------------------------- ASDF images

TEST(ASDFimage, IsAFITSimgAndAnASDFcube)
{
    TempFile path(".asdf");
    writeRoman(path.path());
    ASDFmanager file(path.path());
    auto image = file.GetImage("roman.data");

    EXPECT_TRUE(std::dynamic_pointer_cast<FITSimg<float>>(image));
    auto asdf = asASDF(image);
    ASSERT_TRUE(asdf);
    EXPECT_EQ(&asdf->Cube(), image.get());
    EXPECT_EQ(asdf->ImagePath(), "roman/data");
    EXPECT_EQ(asdf->CompanionNames(), (std::vector<std::string>{"err", "dq", "amp33"}));
    for(const char* name : {"err", "dq", "amp33"})
        EXPECT_FALSE(asdf->Companion(name)->IsLoaded()) << name;
    EXPECT_EQ(asdf->Companion("nothing"), nullptr);

    EXPECT_EQ(image->Size(1), kNX);
    EXPECT_EQ(image->Size(2), kNY);
    EXPECT_DOUBLE_EQ(image->at(5 * kNX + 7), dataAt(7, 5));

    // The FITS side works as for any image
    TempFile fits(".fits");
    image->Write(fits.path());
    FITSmanager back(fits.path());
    EXPECT_DOUBLE_EQ(back.GetPrimary()->at(5 * kNX + 7), dataAt(7, 5));
}

TEST(ASDFimage, WindowCutsTheCompanionsWithoutReadingThem)
{
    TempFile in(".asdf"), out(".asdf");
    writeRoman(in.path());
    ASDFmanager file(in.path());
    auto image = file.GetImage("roman.data");

    auto cut = image->Window(10, 20, 30, 16);
    auto asdf = asASDF(cut);
    ASSERT_TRUE(asdf);
    EXPECT_EQ(cut->Size(1), 30u);
    EXPECT_EQ(cut->Size(2), 16u);
    EXPECT_DOUBLE_EQ(cut->at(0), dataAt(10, 20));

    auto err = asdf->Companion("err");
    ASSERT_TRUE(err);
    EXPECT_TRUE(err->HasView());
    EXPECT_FALSE(err->IsLoaded());
    EXPECT_EQ(err->Axes(), (std::vector<size_t>{30, 16}));
    EXPECT_FLOAT_EQ(static_cast<float>(err->Load()->at(1 * 30 + 2)), errAt(12, 21));
    EXPECT_EQ(asdf->Companion("amp33")->Axes(), (std::vector<size_t>{4, kNY}));
    EXPECT_EQ(asdf->Meta().FindPath("roman/meta/wcs"), nullptr);
    ASSERT_EQ(asdf->PendingHistory().size(), 1u);
    EXPECT_NE(asdf->PendingHistory()[0].find("Window"), std::string::npos);

    // The image it was cut from is unchanged
    EXPECT_FALSE(asASDF(image)->Companion("err")->HasView());
    EXPECT_NE(asASDF(image)->Meta().FindPath("roman/meta/wcs"), nullptr);

    asdf->WriteASDF(out.path());
    ASDFmanager back(out.path());
    ASDFnode tree = back.GetTree();
    for(const char* name : {"data", "err", "dq"})
        EXPECT_EQ(tree.At(std::string("roman/") + name).Array()->Axes(), (std::vector<size_t>{30, 16})) << name;
    EXPECT_EQ(tree.At("roman/amp33").Array()->Axes(), (std::vector<size_t>{4, kNY}));
    EXPECT_DOUBLE_EQ(tree.At("roman/dq").Array()->ValueAt({3, 2}), dqAt(13, 22));
    EXPECT_DOUBLE_EQ(back.GetImage("roman.data")->at(0), dataAt(10, 20));
    EXPECT_EQ(tree.FindPath("roman/meta/wcs"), nullptr);
    EXPECT_EQ(tree.At("roman/meta/exposure/start_time").Tag(), "tag:stsci.edu:asdf/time/time-1.4.0");

    const ASDFnode& entries = tree.At("history/entries");
    ASSERT_EQ(entries.Size(), 2u);
    EXPECT_NE(entries[1].At("description").Get<std::string>().find("Window"), std::string::npos);
    EXPECT_EQ(entries[1].At("software/name").Get<std::string>(), "asdf2fits");
    EXPECT_FALSE(entries[1].At("time").Get<std::string>().empty());
}

TEST(ASDFimage, ResizeCutsTheCompanionsToo)
{
    TempFile in(".asdf");
    writeRoman(in.path());
    ASDFmanager file(in.path());
    auto image = file.GetImage("roman.data");
    image->Resize(4, 6, 32, 24);

    auto asdf = asASDF(image);
    EXPECT_EQ(asdf->Companion("dq")->Axes(), (std::vector<size_t>{32, 24}));
    EXPECT_DOUBLE_EQ(asdf->Companion("dq")->ValueAt({0, 0}), dqAt(4, 6));
    EXPECT_DOUBLE_EQ(image->at(0), dataAt(4, 6));
}

TEST(ASDFimage, LayerTakesTheSameLayerOfTheCompanions)
{
    TempFile in(".asdf"), out(".asdf");
    std::vector<int16_t> cube;
    std::vector<double> var;
    for(int z = 0; z < 3; z++)
        for(int y = 0; y < 5; y++)
            for(int x = 0; x < 4; x++)
            {
                cube.push_back(static_cast<int16_t>(z * 100 + y * 10 + x));
                var.push_back((z * 100 + y * 10 + x) * 0.5);
            }
    writeTestFile(in.path(),
        "sci:\n"
        "  data: !core/ndarray-1.1.0 {source: 0, datatype: int16, byteorder: little, shape: [3, 5, 4]}\n"
        "  var: !core/ndarray-1.1.0 {source: 1, datatype: float64, byteorder: little, shape: [3, 5, 4]}\n",
        {blockOf(cube), blockOf(var, "zlib")});

    ASDFmanager file(in.path());
    auto layer = file.GetImage("sci.data")->Layer(1);
    auto asdf = asASDF(layer);
    ASSERT_TRUE(asdf);
    EXPECT_EQ(layer->GetDimension(), 2u);
    EXPECT_DOUBLE_EQ(layer->at(0), 100.);

    auto companion = asdf->Companion("var");
    EXPECT_EQ(companion->Axes(), (std::vector<size_t>{4, 5}));
    EXPECT_FALSE(companion->IsLoaded());
    EXPECT_DOUBLE_EQ(companion->Load()->at(13), 131 * 0.5);   // x = 1, y = 3

    asdf->WriteASDF(out.path());
    ASDFmanager back(out.path());
    EXPECT_EQ(back.GetArray("sci/var")->Shape(), (std::vector<uint64_t>{5, 4}));
    EXPECT_DOUBLE_EQ(back.GetArray("sci/var")->ValueAt({1, 3}), 131 * 0.5);
}

TEST(ASDFimage, RebinDropsTheCompanionsOfTheImageShape)
{
    TempFile in(".asdf");
    writeRoman(in.path());
    ASDFmanager file(in.path());
    auto binned = file.GetImage("roman.data")->Rebin({2, 2});
    auto asdf = asASDF(binned);
    ASSERT_TRUE(asdf);
    EXPECT_EQ(binned->Size(1), kNX / 2);
    EXPECT_EQ(asdf->CompanionNames(), (std::vector<std::string>{"amp33"}));
    EXPECT_EQ(asdf->Meta().FindPath("roman/err"), nullptr);
    EXPECT_EQ(asdf->Meta().FindPath("roman/meta/wcs"), nullptr);
    ASSERT_EQ(asdf->PendingHistory().size(), 1u);
    EXPECT_NE(asdf->PendingHistory()[0].find("dropped: err, dq"), std::string::npos);

    // FITSimg gives the result a header of its own: the metadata must survive anyway
    EXPECT_EQ(binned->HDU().GetValueForKey("META.INSTRUMENT.NAME"), "WFI");
    TempFile out(".asdf");
    asdf->WriteASDF(out.path());
    ASDFnode tree = ASDFmanager(out.path()).GetTree();
    EXPECT_EQ(tree.At("roman/meta/instrument/name").Get<std::string>(), "WFI");
    EXPECT_EQ(tree.At("roman/meta/exposure/nresultants").Get<int>(), 10);
    EXPECT_EQ(tree.At("roman/data").Array()->Axes(), (std::vector<size_t>{kNX / 2, kNY / 2}));
}

TEST(ASDFimage, ShapeChangedOutsideTheHooksDropsTheCompanions)
{
    TempFile in(".asdf"), out(".asdf");
    writeRoman(in.path());
    ASDFmanager file(in.path());
    auto image = file.GetImage("roman.data");

    // A call through the FITSimg type, which ASDFimg can't see
    FITSimg<float>& plain = dynamic_cast<FITSimg<float>&>(*image);
    plain.FITSimg<float>::AddLayer(static_cast<const FITSimg<float>&>(*ramp<float>({kNX, kNY})));
    ASSERT_EQ(image->GetDimension(), 3u);

    auto asdf = asASDF(image);
    EXPECT_EQ(asdf->CompanionNames(), (std::vector<std::string>{"amp33"}));
    asdf->WriteASDF(out.path());
    ASDFnode tree = ASDFmanager(out.path()).GetTree();
    EXPECT_EQ(tree.FindPath("roman/err"), nullptr);
    EXPECT_NE(tree.FindPath("roman/amp33"), nullptr);
    EXPECT_EQ(tree.FindPath("roman/meta/wcs"), nullptr);
    EXPECT_NE(tree.At("history/entries/1/description").Get<std::string>().find("reshaped"), std::string::npos);
}

TEST(ASDFimage, MaskFromDQReadsTheFlagsAndFreesThem)
{
    TempFile in(".asdf");
    writeRoman(in.path());
    ASDFmanager file(in.path());
    auto image = file.GetImage("roman.data");
    auto asdf = asASDF(image);

    EXPECT_EQ(asdf->MaskFromDQ(1), 9 * kNY);   // x % 7 == 3: 9 pixels per row
    EXPECT_TRUE(image->Masked(size_t(5 * kNX + 3)));
    EXPECT_FALSE(image->Masked(size_t(0)));
    EXPECT_FALSE(asdf->Companion("dq")->IsLoaded());

    EXPECT_EQ(asdf->MaskFromDQ(4), kNX);       // all of row 5
    EXPECT_TRUE(image->Masked(size_t(5 * kNX)));
    EXPECT_THROW(asdf->MaskFromDQ(1, "err"), FITSexception);
    EXPECT_THROW(asdf->MaskFromDQ(1, "amp33"), FITSexception);
    EXPECT_THROW(asdf->MaskFromDQ(1, "nothing"), FITSexception);
}

TEST(ASDFimage, AsdfLikeGivesAResultTheTree)
{
    TempFile in(".asdf");
    writeRoman(in.path());
    ASDFmanager file(in.path());
    auto model = asASDF(file.GetImage("roman.data"));

    // Kept elsewhere: copied
    auto stack = ramp<float>({kNX, kNY});
    auto like = asASDF(asdfLike(stack, *model, "stacked"));
    ASSERT_TRUE(like);
    EXPECT_DOUBLE_EQ(like->Cube().at(3), stack->at(3));
    EXPECT_EQ(like->CompanionNames(), (std::vector<std::string>{"amp33"}));
    EXPECT_NE(like->Meta().FindPath("roman/meta/wcs"), nullptr);
    EXPECT_EQ(like->PendingHistory().back(), "stacked; companions dropped: err, dq");

    // Of another shape, and given away: moved, without its GWCS
    auto smaller = asASDF(asdfLike(ramp<float>({8, 8}), *model, "binned"));
    ASSERT_TRUE(smaller);
    EXPECT_EQ(smaller->Cube().Size(1), 8u);
    EXPECT_EQ(smaller->Meta().FindPath("roman/meta/wcs"), nullptr);

    // A result made with a header of its own keeps the metadata, in the header and in the file
    EXPECT_EQ(like->Cube().HDU().GetValueForKey("META.INSTRUMENT.NAME"), "WFI");
    TempFile out(".asdf");
    like->WriteASDF(out.path());
    ASDFnode tree = ASDFmanager(out.path()).GetTree();
    EXPECT_EQ(tree.At("roman/meta/instrument/name").Get<std::string>(), "WFI");
    EXPECT_EQ(tree.At("roman/meta/exposure/start_time").Tag(), "tag:stsci.edu:asdf/time/time-1.4.0");
    EXPECT_NE(tree.FindPath("roman/amp33"), nullptr);
    EXPECT_EQ(tree.FindPath("roman/err"), nullptr);
}

TEST(ASDFmanager, GetImageWindowReadsOnlyTheWindow)
{
    TempFile in(".asdf");
    writeRoman(in.path());
    ASDFmanager file(in.path());

    auto cut = file.GetImageWindow("roman.data", {8, 4}, {16, 8});
    EXPECT_EQ(cut->Size(1), 16u);
    EXPECT_EQ(cut->Size(2), 8u);
    EXPECT_DOUBLE_EQ(cut->at(0), dataAt(8, 4));
    EXPECT_DOUBLE_EQ(cut->at(15 + 7 * 16), dataAt(23, 11));
    EXPECT_EQ(cut->HDU().GetInt64ValueForKey("NAXIS1"), 16);
    EXPECT_EQ(cut->HDU().GetInt64ValueForKey("NAXIS2"), 8);

    auto asdf = asASDF(cut);
    EXPECT_EQ(asdf->Companion("err")->Axes(), (std::vector<size_t>{16, 8}));
    EXPECT_EQ(asdf->Meta().FindPath("roman/meta/wcs"), nullptr);
    EXPECT_NE(asdf->PendingHistory().at(0).find("GetImageWindow"), std::string::npos);

    // To the end of the axes
    EXPECT_EQ(file.GetImageWindow("roman.data", {60})->Size(1), 4u);
    EXPECT_THROW(file.GetImageWindow("roman.data", {kNX, 0}), FITSexception);
}

TEST(ASDFmanager, GetImageWindowShiftsTheFITSWCS)
{
    TempFile in(".asdf");
    std::vector<float> data(10 * 8, 1.f);
    writeTestFile(in.path(),
        "sci:\n"
        "  CTYPE1: RA---TAN\n  CTYPE2: DEC--TAN\n  CRPIX1: 10.0\n  CRPIX2: 20.0\n"
        "  CRVAL1: 10.0\n  CRVAL2: 20.0\n  CDELT1: -0.001\n  CDELT2: 0.001\n"
        "  data: !core/ndarray-1.1.0 {source: 0, datatype: float32, byteorder: little, shape: [8, 10]}\n",
        {blockOf(data)});
    ASDFmanager file(in.path());
    auto cut = file.GetImageWindow("sci.data", {3, 5}, {4, 2});
    EXPECT_DOUBLE_EQ(cut->HDU().GetDoubleValueForKey("CRPIX1"), 7.);
    EXPECT_DOUBLE_EQ(cut->HDU().GetDoubleValueForKey("CRPIX2"), 15.);
}

// ---------------------------------------------------------------------- Header and tree

TEST(ASDFsync, HeaderChangesReachTheTree)
{
    TempFile in(".asdf"), out(".asdf");
    writeRoman(in.path());
    ASDFmanager file(in.path());
    auto image = file.GetImage("roman.data");

    FITShdu& hdu = image->HDU();
    hdu.ValueForKey("META.EXPOSURE.EXPOSURE_TIME", 200.5);
    hdu.ValueForKey("META.EXPOSURE.START_TIME", std::string("2028-01-01T00:00:00.000"), fChar);
    hdu.DeleteKey("META.INSTRUMENT.DETECTOR");
    hdu.ValueForKey("OBSERVER", std::string("me"), fChar);
    hdu.ValueForKey("META.EXPOSURE.READS", int32_t(7));
    asASDF(image)->WriteASDF(out.path());

    ASDFmanager back(out.path());
    ASDFnode tree = back.GetTree();
    EXPECT_DOUBLE_EQ(tree.At("roman/meta/exposure/exposure_time").Get<double>(), 200.5);
    EXPECT_EQ(tree.At("roman/meta/exposure/start_time").Get<std::string>(), "2028-01-01T00:00:00.000");
    EXPECT_EQ(tree.At("roman/meta/exposure/start_time").Tag(), "tag:stsci.edu:asdf/time/time-1.4.0");
    EXPECT_EQ(tree.FindPath("roman/meta/instrument/detector"), nullptr);
    EXPECT_EQ(tree.At("roman/meta/instrument/name").Get<std::string>(), "WFI");
    EXPECT_EQ(tree.At("roman/OBSERVER").Get<std::string>(), "me");
    EXPECT_EQ(tree.At("roman/meta/exposure/reads").Get<int>(), 7);
    for(const char* key : {"BITPIX", "NAXIS", "NAXIS1", "BSCALE", "BZERO", "EXTNAME"})
        EXPECT_FALSE(tree.At("roman").Contains(key)) << key;

    // Read back, the header is the one that was written
    auto header = back.GetHeader("roman.data");
    EXPECT_EQ(header->GetValueForKey("OBSERVER"), "me");
    EXPECT_EQ(header->GetInt64ValueForKey("META.EXPOSURE.READS"), 7);
    EXPECT_FALSE(header->Exists("META.INSTRUMENT.DETECTOR"));
}

TEST(ASDFsync, SetAndRefreshKeepTheHeaderInStep)
{
    TempFile in(".asdf");
    writeRoman(in.path());
    ASDFmanager file(in.path());
    auto image = file.GetImage("roman.data");
    auto asdf = asASDF(image);

    asdf->Set("roman/meta/instrument/name", ASDFnode("NISP"));
    EXPECT_EQ(image->HDU().GetValueForKey("META.INSTRUMENT.NAME"), "NISP");
    EXPECT_TRUE(asdf->Erase("roman/meta/exposure/nresultants"));
    EXPECT_FALSE(image->HDU().Exists("META.EXPOSURE.NRESULTANTS"));

    // A change of the header isn't lost when the tree is read again
    image->HDU().ValueForKey("OBSERVER", std::string("me"), fChar);
    asdf->Meta().SetPath("roman/meta/new", 1.5);
    asdf->RefreshHDU();
    EXPECT_DOUBLE_EQ(image->HDU().GetDoubleValueForKey("META.NEW"), 1.5);
    EXPECT_EQ(asdf->Meta().At("roman/OBSERVER").Get<std::string>(), "me");
    EXPECT_EQ(image->HDU().GetValueForKey("OBSERVER"), "me");
}

// ---------------------------------------------------------------------- Writing

TEST(ASDFwrite, WhatIsntChangedComesBackAsItWas)
{
    TempFile in(".asdf"), out(".asdf");
    writeRoman(in.path());
    {
        ASDFmanager file(in.path());
        asASDF(file.GetImage("roman.data"))->WriteASDF(out.path());
    }

    ASDFmanager a(in.path()), b(out.path());
    expectSameTree(a.GetTree().At("roman"), b.GetTree().At("roman"));

    const ASDFnode entries = b.GetTree().At("history/entries");
    ASSERT_EQ(entries.Size(), 1u);
    EXPECT_EQ(entries[0].At("description").Get<std::string>(), "made by the test");
    EXPECT_EQ(entries[0].At("time").Get<std::string>(), "2020-01-01T00:00:00");

    // Arrays that weren't loaded are copied as they were, still compressed
    for(const char* name : {"err", "dq", "amp33"})
    {
        const auto before = rawBlock(in.path(), std::string("roman/") + name);
        const auto after  = rawBlock(out.path(), std::string("roman/") + name);
        EXPECT_FALSE(before.first.empty()) << name;
        EXPECT_EQ(before.first, after.first) << name;
        EXPECT_TRUE(before.second == after.second) << name;
    }
}

TEST(ASDFwrite, RomanFixturesComeBackAsTheyWere)
{
    for(const char* name : {"roman_l2_wcs", "roman_wcs"})
    {
        TempFile out(".asdf");
        ASDFmanager file(fixture(name));
        const ASDFnode tree = file.GetTree();
        WriteASDF(tree, out.path());

        ASDFnode back = ASDFmanager(out.path()).GetTree();
        for(size_t i = 0; i < tree.Size(); i++)
            if(tree.Keys()[i] != "asdf_library" && tree.Keys()[i] != "history")
                expectSameTree(tree.Items()[i], back.At(tree.Keys()[i]), std::string(name) + "/" + tree.Keys()[i]);
    }
}

TEST(ASDFwrite, CompressionOptions)
{
    TempFile in(".asdf");
    writeRoman(in.path());
    ASDFmanager file(in.path());
    auto image = file.GetImage("roman.data");
    auto asdf = asASDF(image);
    asdf->Companion("err")->Load();

    auto written = [&](const std::string& compression)
    {
        TempFile out(".asdf");
        ASDFwriteOptions options;
        options.compression = compression;
        asdf->WriteASDF(out.path(), options);
        ASDFnode tree = ASDFmanager(out.path()).GetTree();
        std::vector<std::string> found;
        for(const char* name : {"data", "err", "dq", "amp33"})
            found.push_back(tree.At(std::string("roman/") + name).Array()->Compression());
        EXPECT_DOUBLE_EQ(tree.At("roman/data").Array()->ValueAt({7, 5}), dataAt(7, 5));
        EXPECT_FLOAT_EQ(static_cast<float>(tree.At("roman/err").Array()->ValueAt({7, 5})), errAt(7, 5));
        return found;
    };

    // Written from memory: data and the loaded err; copied: dq and amp33
    EXPECT_EQ(written("zlib"), (std::vector<std::string>{"zlib", "zlib", "zlib", "bzp2"}));
    EXPECT_EQ(written("keep"), (std::vector<std::string>{"", "lz4", "zlib", "bzp2"}));
    EXPECT_EQ(written(""),     (std::vector<std::string>{"", "", "zlib", "bzp2"}));
    EXPECT_EQ(written("lz4"),  (std::vector<std::string>{"lz4", "lz4", "zlib", "bzp2"}));

    TempFile bad(".asdf");
    ASDFwriteOptions options;
    options.compression = "gzip";
    EXPECT_THROW(asdf->WriteASDF(bad.path(), options), FITSexception);
}

TEST(ASDFwrite, OverTheFileItWasReadFrom)
{
    TempFile in(".asdf");
    writeRoman(in.path());
    std::shared_ptr<FITScube> image;
    {
        ASDFmanager file(in.path());
        image = file.GetImage("roman.data");
    }
    image->SetPixelValue(-1.f, size_t(0));
    EXPECT_THROW(asASDF(image)->WriteASDF(in.path()), FITSexception);

    ASDFwriteOptions options;
    options.replace = true;
    asASDF(image)->WriteASDF(in.path(), options);

    ASDFmanager back(in.path());
    EXPECT_DOUBLE_EQ(back.GetImage("roman.data")->at(0), -1.);
    EXPECT_DOUBLE_EQ(back.GetImage("roman.data")->at(1), dataAt(1, 0));
    EXPECT_FLOAT_EQ(static_cast<float>(back.GetArray("roman/err")->ValueAt({2, 3})), errAt(2, 3));
    EXPECT_DOUBLE_EQ(back.GetArray("roman/dq")->ValueAt({3, 5}), dqAt(3, 5));
}

TEST(ASDFwrite, HistoryEntriesAreAdded)
{
    TempFile in(".asdf"), out(".asdf");
    writeRoman(in.path());
    ASDFmanager file(in.path());
    auto asdf = asASDF(file.GetImage("roman.data"));
    asdf->AddHistory("calibrated");
    ASDFwriteOptions options;
    options.history = "saved by the test";
    asdf->WriteASDF(out.path(), options);

    const ASDFnode entries = ASDFmanager(out.path()).GetTree().At("history/entries");
    ASSERT_EQ(entries.Size(), 3u);
    EXPECT_EQ(entries[0].At("time").Get<std::string>(), "2020-01-01T00:00:00");
    EXPECT_EQ(entries[1].At("description").Get<std::string>(), "calibrated");
    EXPECT_EQ(entries[2].At("description").Get<std::string>(), "saved by the test");
    EXPECT_EQ(entries[2].Tag(), "tag:stsci.edu:asdf/core/history_entry-1.0.0");
    EXPECT_EQ(entries[2].At("software/name").Get<std::string>(), "asdf2fits");
    EXPECT_EQ(entries[2].At("time").Get<std::string>().size(), 19u);
}

TEST(ASDFwrite, TreeBuiltFromScratch)
{
    TempFile out(".asdf");
    ASDFnode tree;
    tree["meta"]["object"] = "M31";
    tree["meta"]["exptime"] = 30.;
    tree["sci"]["data"] = ASDFnode(ASDFarray::FromCube(ramp<float>({6, 4})));
    WriteASDF(tree, out.path());

    ASDFmanager back(out.path());
    auto image = back.GetImage("sci.data");
    EXPECT_DOUBLE_EQ(image->at(5), 6.);
    EXPECT_EQ(image->HDU().GetValueForKey("META.OBJECT"), "M31");
    EXPECT_DOUBLE_EQ(image->HDU().GetDoubleValueForKey("META.EXPTIME"), 30.);
}

TEST(ASDFwrite, AstropyLayoutKeepsItsCards)
{
    TempFile in(".asdf"), out(".asdf");
    auto source = ramp<float>({6, 4});
    source->HDU().ValueForKey("OBSERVER", std::string("someone"), fChar, "who");
    fits2asdf(*source, in.path());

    ASDFmanager file(in.path());
    auto image = file.GetPrimary();
    ASSERT_TRUE(asASDF(image));
    EXPECT_TRUE(asASDF(image)->CompanionNames().empty());
    image->HDU().ValueForKey("OBSERVER", std::string("me"), fChar);
    image->HDU().ValueForKey("NEWKEY", int32_t(3));
    asASDF(image->Window(1, 1, 3, 2))->WriteASDF(out.path());

    ASDFmanager back(out.path());
    auto header = back.GetPrimaryHeader();
    EXPECT_EQ(header->GetValueForKey("OBSERVER"), "me");
    EXPECT_EQ(header->GetInt64ValueForKey("NEWKEY"), 3);
    EXPECT_EQ(header->GetInt64ValueForKey("NAXIS1"), 3);
    EXPECT_EQ(header->GetInt64ValueForKey("NAXIS2"), 2);
    EXPECT_DOUBLE_EQ(back.GetPrimary()->at(0), source->at(1 * 6 + 1));

    ASDFnode tree = back.GetTree();
    EXPECT_TRUE(tree.At("fits").HasTag("/fits/fits-"));
    EXPECT_EQ(tree.At("fits/0/header/0/0").Get<std::string>(), "SIMPLE");
}

TEST(ASDFwrite, Fits2asdfWritesAnASDFImageInItsOwnLayout)
{
    TempFile in(".asdf"), out(".asdf");
    writeRoman(in.path());
    ASDFmanager file(in.path());
    fits2asdf(*file.GetImage("roman.data"), out.path());
    ASDFnode tree = ASDFmanager(out.path()).GetTree();
    EXPECT_EQ(tree.FindPath("fits"), nullptr);
    EXPECT_EQ(tree.At("roman/meta/instrument/name").Get<std::string>(), "WFI");
}

// ---------------------------------------------------------------------- Array handles

TEST(ASDFarray, LoadReadReleaseSetAndClone)
{
    TempFile in(".asdf");
    writeRoman(in.path());
    ASDFmanager file(in.path());
    auto err = file.GetArray("roman/err");
    EXPECT_FALSE(err->IsLoaded());
    EXPECT_THROW(file.GetArray("roman/meta"), FITSexception);
    EXPECT_THROW(file.GetArray("roman/none"), FITSexception);

    auto pixels = err->Load();
    EXPECT_TRUE(err->IsLoaded());
    EXPECT_EQ(err->Load(), pixels);

    auto tile = err->Read({1, 2}, {3, 4});
    EXPECT_EQ(tile->Size(1), 3u);
    EXPECT_EQ(tile->Size(2), 4u);
    EXPECT_FLOAT_EQ(static_cast<float>(tile->at(0)), errAt(1, 2));
    EXPECT_THROW(err->Read({kNX}), FITSexception);

    auto copy = err->Clone();
    EXPECT_TRUE(copy->IsLoaded());
    EXPECT_NE(copy->Load(), pixels);
    EXPECT_TRUE(err->Release());
    EXPECT_FALSE(err->IsLoaded());

    auto window = err->Window({2, 3}, {4, 5})->Window({1, 1}, {2, 2});
    EXPECT_EQ(window->Axes(), (std::vector<size_t>{2, 2}));
    EXPECT_FLOAT_EQ(static_cast<float>(window->ValueAt({1, 1})), errAt(4, 5));
    EXPECT_FALSE(window->IsLoaded());

    // The array outlives its manager
    file.Close();
    EXPECT_FLOAT_EQ(static_cast<float>(window->Load()->at(0)), errAt(3, 4));

    auto memory = ASDFarray::FromCube(ramp<int16_t>({5, 4}));
    EXPECT_FALSE(memory->IsInFile());
    EXPECT_EQ(memory->DataType(), "int16");
    EXPECT_EQ(memory->Axes(), (std::vector<size_t>{5, 4}));
    EXPECT_EQ(memory->Bytes(), 40u);
    EXPECT_FALSE(memory->Release());
    EXPECT_TRUE(memory->IsLoaded());
}

// libasdf writes over a file without truncating it: a longer file replaced by a shorter one kept its end,
// and its block index, behind.
TEST(Write, ReplacingALongerFileLeavesNothingOfIt)
{
    TempFile path(".asdf"), fresh(".asdf");
    fits2asdf(*ramp<double>({200, 200}), path.path());
    auto small = ramp<float>({4, 3});
    fits2asdf(*small, path.path(), true);
    fits2asdf(*small, fresh.path());

    EXPECT_EQ(std::filesystem::file_size(path.path()), std::filesystem::file_size(fresh.path()));
    ASDFmanager back(path.path());
    EXPECT_EQ(back.GetPrimary()->Nelements(), 12u);
}

TEST(ASDFmanager, LargeCompressedArraysCanBeDecompressedToDisk)
{
    TempFile in(".asdf");
    writeRoman(in.path());
    ASDFreadOptions options;
    options.maxMemoryBytes = 1024;      // smaller than every array: all go through a temporary file
    options.tmpDir = ASDF2FITS_TEST_DIR;    // created by the TempFile above
    ASDFmanager file(in.path(), options);

    auto dq = file.GetArray("roman/dq");
    EXPECT_DOUBLE_EQ(dq->ValueAt({3, 5}), dqAt(3, 5));
    auto err = file.GetArray("roman/err")->Load();
    EXPECT_FLOAT_EQ(static_cast<float>(err->at(7 * kNX + 9)), errAt(9, 7));
    EXPECT_DOUBLE_EQ(file.GetArray("roman/amp33")->ValueAt({2, 40}), amp33At(2, 40));
}

// ---------------------------------------------------------------------- Multi-HDU file, libasdf on the other side
//
// The file of writeMultiHDU() is made with the libasdf API alone: one HDU per pixel type, of 1 to 4 axes, little
// or big endian, compressed or not, plus a mapping whose image has companions. ASDFcube must read exactly what
// was put in it, and what ASDFcube writes must read back, with libasdf alone, as expected.

template<typename T> class LibasdfMultiHDU : public testing::Test {};
TYPED_TEST_SUITE(LibasdfMultiHDU, PixelTypes);

TYPED_TEST(LibasdfMultiHDU, ASDFcubeReadsWhatLibasdfWrote)
{
    using T = TypeParam;
    const ArraySpec spec = specOf<T>();
    const std::string name = spec.group + ".data";
    TempFile path(".asdf");
    ASSERT_NO_FATAL_FAILURE(writeMultiHDU(path.path()));

    ASDFmanager file(path.path());
    ASSERT_EQ(file.NumberOfHeader(), kMultiHDUs);

    auto cube = file.GetImage(name);
    auto image = std::dynamic_pointer_cast<ASDFimg<T>>(cube);
    ASSERT_NE(image, nullptr) << name << " isn't an ASDFimg of its type";
    EXPECT_EQ(image->ImagePath(), spec.group + "/data");
    EXPECT_TRUE(image->CompanionNames().empty());

    // Axes: NAXIS1 is the last length of the shape
    const std::vector<size_t> axes = axesOf(spec.shape);
    ASSERT_EQ(cube->GetDimension(), axes.size());
    for(size_t k = 0; k < axes.size(); k++)
        EXPECT_EQ(cube->Size(k + 1), axes[k]) << "NAXIS" << k + 1;

    // Pixels: exactly the values written, whatever the byte order and the compression
    const std::vector<T> expected = pixelsOf<T>(spec);
    const std::valarray<T>* pixels = cube->template GetData<T>();
    ASSERT_NE(pixels, nullptr);
    ASSERT_EQ(pixels->size(), expected.size());
    for(size_t i = 0; i < expected.size(); i++)
        ASSERT_EQ((*pixels)[i], expected[i]) << name << " [" << i << "]";
    EXPECT_EQ(cube->GetBitPerPixel(), FITSimg<T>(std::vector<size_t>{1}).GetBitPerPixel());

    // Header: structure, and the scalars of the mapping of the image
    auto header = file.GetHeader(name);
    EXPECT_EQ(header->GetValueForKey("EXTNAME"), name);
    EXPECT_EQ(header->GetInt32ValueForKey("NAXIS"), static_cast<int32_t>(axes.size()));
    for(size_t k = 0; k < axes.size(); k++)
        EXPECT_EQ(header->GetUInt64ValueForKey("NAXIS" + std::to_string(k + 1)), axes[k]);
    EXPECT_EQ(header->GetValueForKey("LABEL"), spec.group + " image");
    EXPECT_EQ(header->GetInt64ValueForKey("RANK"), static_cast<int64_t>(axes.size()));
    EXPECT_EQ(cube->HDU().GetValueForKey("LABEL"), spec.group + " image");

    // Tree: the array as described in the file
    const ASDFnode tree = file.GetTree();
    auto array = tree.At(spec.group + "/data").Array();
    ASSERT_TRUE(array);
    EXPECT_EQ(array->DataType(), asdfNameOf<T>());
    EXPECT_EQ(array->Shape(), spec.shape);
    EXPECT_EQ(array->Axes(), axes);
    EXPECT_EQ(array->Compression(), spec.compression);
}

TEST(LibasdfMultiHDU, ASDFcubeReadsTheFileLevelContent)
{
    TempFile path(".asdf");
    ASSERT_NO_FATAL_FAILURE(writeMultiHDU(path.path()));
    ASDFmanager file(path.path());
    ASSERT_EQ(file.NumberOfHeader(), kMultiHDUs);

    // Every ndarray is an HDU, named after its path
    std::vector<std::string> names;
    for(int i = 1; i <= file.NumberOfHeader(); i++)
        names.push_back(file.GetHeaderAtIndex(i)->GetValueForKey("EXTNAME"));
    std::vector<std::string> expected = {"sci.data", "sci.err", "sci.dq"};
    forEachPixelType([&]<typename T>() { expected.push_back(specOf<T>().group + ".data"); });
    std::sort(names.begin(), names.end());
    std::sort(expected.begin(), expected.end());
    EXPECT_EQ(names, expected);

    // The scalars of the root are keywords of the primary HDU
    auto primary = file.GetPrimaryHeader();
    EXPECT_EQ(primary->GetValueForKey("OBJECT"), "NGC 1300");
    EXPECT_DOUBLE_EQ(primary->GetDoubleValueForKey("EXPTIME"), 87.5);
    EXPECT_EQ(primary->GetInt64ValueForKey("NCOMBINE"), -12);
    EXPECT_TRUE(primary->GetBoolValueForKey("CALIBRATED"));

    // The image of sci has its siblings as companions, read only when asked
    auto sci = asASDF(file.GetImage("sci.data"));
    ASSERT_TRUE(sci);
    EXPECT_EQ(sci->CompanionNames(), (std::vector<std::string>{"err", "dq"}));
    EXPECT_EQ(sci->Cube().HDU().GetValueForKey("TELESCOPE"), "Euclid");
    EXPECT_EQ(sci->Cube().Size(1), kSciNX);
    EXPECT_EQ(sci->Cube().Size(2), kSciNY);
    EXPECT_DOUBLE_EQ(sci->Cube().at(3 * kSciNX + 5), sciDataAt(5, 3));

    auto err = sci->Companion("err");
    ASSERT_TRUE(err);
    EXPECT_FALSE(err->IsLoaded());
    EXPECT_EQ(err->Compression(), "lz4");
    EXPECT_EQ(err->DataType(), "float32");
    EXPECT_FLOAT_EQ(static_cast<float>(err->ValueAt({5, 3})), sciErrAt(5, 3));     // big endian, compressed

    EXPECT_EQ(sci->MaskFromDQ(1), kSciNY);     // column x = 2
    EXPECT_EQ(sci->MaskFromDQ(2), kSciNX);     // row y = 4
    EXPECT_TRUE(sci->Cube().Masked(size_t(4 * kSciNX + 7)));
    EXPECT_FALSE(sci->Cube().Masked(size_t(3 * kSciNX + 5)));

    // Companions are HDUs of their own too
    auto dq = file.GetImage("sci.dq");
    ASSERT_NE(dq->GetData<uint16_t>(), nullptr);
    EXPECT_EQ(dq->UShortValueAtPixel(size_t(4 * kSciNX + 2)), 3);
    auto errImage = file.GetImage("sci.err");
    ASSERT_NE(errImage->GetData<float>(), nullptr);
    for(size_t y = 0; y < kSciNY; y++)
        for(size_t x = 0; x < kSciNX; x++)
            ASSERT_EQ(errImage->FloatValueAtPixel(y * kSciNX + x), sciErrAt(x, y)) << x << "," << y;
}

TYPED_TEST(LibasdfMultiHDU, LibasdfReadsWhatASDFcubeWrote)
{
    using T = TypeParam;
    const ArraySpec spec = specOf<T>();
    TempFile in(".asdf"), out(".asdf");
    ASSERT_NO_FATAL_FAILURE(writeMultiHDU(in.path()));

    // Change one pixel in three, and a keyword
    std::vector<T> expected = pixelsOf<T>(spec);
    const std::vector<T> changed = pixelsOf<T>(spec, spec.salt + 50);
    {
        ASDFmanager file(in.path());
        auto image = std::dynamic_pointer_cast<ASDFimg<T>>(file.GetImage(spec.group + ".data"));
        ASSERT_NE(image, nullptr);
        std::valarray<T>& pixels = *image->template GetData<T>();
        for(size_t i = 2; i < pixels.size(); i += 3)
            pixels[i] = expected[i] = changed[i];
        image->HDU().ValueForKey("LABEL", std::string("changed"), fChar);
        image->WriteASDF(out.path());
    }

    // The image, as libasdf reads it: type, shape, compression kept, pixels
    LibasdfArray back;
    ASSERT_NO_FATAL_FAILURE(readNdarray(out.path(), spec.group + "/data", back));
    EXPECT_EQ(back.type, asdfTypeOf<T>());
    EXPECT_EQ(back.shape, spec.shape);
    EXPECT_EQ(back.compression, spec.compression);
    EXPECT_TRUE(back.values<T>() == expected);

    // The arrays of the other HDUs, still compressed, come back as they were
    for(const char* other : {"sci/err", "sci/dq"})
    {
        const auto before = rawBlock(in.path(), other);
        const auto after  = rawBlock(out.path(), other);
        EXPECT_EQ(before.first, after.first) << other;
        EXPECT_TRUE(before.second == after.second) << other;
    }
    LibasdfArray err;
    ASSERT_NO_FATAL_FAILURE(readNdarray(out.path(), "sci/err", err));
    EXPECT_EQ(err.type, ASDF_DATATYPE_FLOAT32);
    EXPECT_FLOAT_EQ(err.values<float>().at(3 * kSciNX + 5), sciErrAt(5, 3));

    // The scalars: the keyword changed, the others as they were
    asdf_file_t* file = asdf_open_file(out.path().c_str(), "r");
    ASSERT_NE(file, nullptr);
    const char* label = nullptr;
    const char* object = nullptr;
    const char* telescope = nullptr;
    double exptime = 0.;
    int64_t ncombine = 0, rank = 0;
    bool calibrated = false;
    EXPECT_EQ(asdf_get_string0(file, (spec.group + "/label").c_str(), &label), ASDF_VALUE_OK);
    EXPECT_STREQ(label, "changed");
    EXPECT_EQ(asdf_get_int64(file, (spec.group + "/rank").c_str(), &rank), ASDF_VALUE_OK);
    EXPECT_EQ(rank, static_cast<int64_t>(spec.shape.size()));
    EXPECT_EQ(asdf_get_string0(file, "object", &object), ASDF_VALUE_OK);
    EXPECT_STREQ(object, "NGC 1300");
    EXPECT_EQ(asdf_get_double(file, "exptime", &exptime), ASDF_VALUE_OK);
    EXPECT_DOUBLE_EQ(exptime, 87.5);
    EXPECT_EQ(asdf_get_int64(file, "ncombine", &ncombine), ASDF_VALUE_OK);
    EXPECT_EQ(ncombine, -12);
    EXPECT_EQ(asdf_get_bool(file, "calibrated", &calibrated), ASDF_VALUE_OK);
    EXPECT_TRUE(calibrated);
    EXPECT_EQ(asdf_get_string0(file, "sci/telescope", &telescope), ASDF_VALUE_OK);
    EXPECT_STREQ(telescope, "Euclid");
    asdf_close(file);

    // No HDU was lost
    EXPECT_EQ(ASDFmanager(out.path()).NumberOfHeader(), kMultiHDUs);
}

TEST(LibasdfMultiHDU, LibasdfReadsCutImagesWrittenByASDFcube)
{
    TempFile in(".asdf"), window(".asdf"), layer(".asdf");
    ASSERT_NO_FATAL_FAILURE(writeMultiHDU(in.path()));
    const ArraySpec u16 = specOf<uint16_t>();      // shape {7, 6}: NAXIS1 = 6, NAXIS2 = 7
    const ArraySpec i16 = specOf<int16_t>();       // shape {3, 4, 5}
    const std::vector<uint16_t> u16Pixels = pixelsOf<uint16_t>(u16);
    const std::vector<int16_t>  i16Pixels = pixelsOf<int16_t>(i16);
    {
        ASDFmanager file(in.path());
        asASDF(file.GetImage("u16.data")->Window(1, 2, 4, 3))->WriteASDF(window.path());
        asASDF(file.GetImage("i16.data")->Layer(1))->WriteASDF(layer.path());
    }

    LibasdfArray cut;
    ASSERT_NO_FATAL_FAILURE(readNdarray(window.path(), "u16/data", cut));
    EXPECT_EQ(cut.type, ASDF_DATATYPE_UINT16);
    ASSERT_EQ(cut.shape, (std::vector<uint64_t>{3, 4}));
    EXPECT_EQ(cut.compression, "lz4");
    const std::vector<uint16_t> cutPixels = cut.values<uint16_t>();
    ASSERT_EQ(cutPixels.size(), 12u);
    for(size_t y = 0; y < 3; y++)
        for(size_t x = 0; x < 4; x++)
            EXPECT_EQ(cutPixels[y * 4 + x], u16Pixels[(y + 2) * 6 + (x + 1)]) << x << "," << y;

    LibasdfArray plane;
    ASSERT_NO_FATAL_FAILURE(readNdarray(layer.path(), "i16/data", plane));
    EXPECT_EQ(plane.type, ASDF_DATATYPE_INT16);
    ASSERT_EQ(plane.shape, (std::vector<uint64_t>{4, 5}));
    EXPECT_EQ(plane.compression, "");
    const std::vector<int16_t> planePixels = plane.values<int16_t>();
    ASSERT_EQ(planePixels.size(), 20u);
    for(size_t i = 0; i < 20; i++)
        EXPECT_EQ(planePixels[i], i16Pixels[20 + i]) << i;

    // The other HDUs are untouched
    LibasdfArray other;
    ASSERT_NO_FATAL_FAILURE(readNdarray(window.path(), "f64/data", other));
    EXPECT_EQ(other.shape, specOf<double>().shape);
    EXPECT_TRUE(other.values<double>() == pixelsOf<double>(specOf<double>()));
}

TEST(LibasdfMultiHDU, LibasdfReadsATreeOfAllPixelTypesWrittenByASDF2fits)
{
    TempFile out(".asdf");
    ASDFnode tree;
    tree["object"] = "M33";
    forEachPixelType([&]<typename T>()
    {
        const ArraySpec spec = specOf<T>();
        auto cube = std::make_shared<FITSimg<T>>(axesOf(spec.shape));
        std::valarray<T>& pixels = *cube->template GetData<T>();
        const std::vector<T> values = pixelsOf<T>(spec);
        for(size_t i = 0; i < values.size(); i++)
            pixels[i] = values[i];
        tree[spec.group]["data"] = ASDFnode(ASDFarray::FromCube(cube));
        tree[spec.group]["label"] = spec.group + " image";
    });
    ASDFwriteOptions options;
    options.compression = "zlib";
    WriteASDF(tree, out.path(), options);

    forEachPixelType([&]<typename T>()
    {
        const ArraySpec spec = specOf<T>();
        LibasdfArray back;
        ASSERT_NO_FATAL_FAILURE(readNdarray(out.path(), spec.group + "/data", back));
        EXPECT_EQ(back.type, asdfTypeOf<T>()) << spec.group;
        EXPECT_EQ(back.shape, spec.shape) << spec.group;
        EXPECT_EQ(back.compression, "zlib") << spec.group;
        EXPECT_TRUE(back.values<T>() == pixelsOf<T>(spec)) << spec.group;
    });

    asdf_file_t* file = asdf_open_file(out.path().c_str(), "r");
    ASSERT_NE(file, nullptr);
    const char* object = nullptr;
    const char* label = nullptr;
    EXPECT_EQ(asdf_get_string0(file, "object", &object), ASDF_VALUE_OK);
    EXPECT_STREQ(object, "M33");
    EXPECT_EQ(asdf_get_string0(file, "f32/label", &label), ASDF_VALUE_OK);
    EXPECT_STREQ(label, "f32 image");
    asdf_close(file);
}

// ---------------------------------------------------------------------- GWCS

namespace
{
    // A transform as the tree holds it
    ASDFnode transformNode(const std::string& kind, const std::string& version = "1.3.0")
    {
        return ASDFnode::Mapping("tag:stsci.edu:asdf/transform/" + kind + "-" + version);
    }

    ASDFnode numberList(std::initializer_list<double> values)
    {
        ASDFnode list = ASDFnode::Sequence();
        for(double value : values)
            list.Append(value);
        return list;
    }

    ASDFnode nodeList(std::initializer_list<ASDFnode> items)
    {
        ASDFnode list = ASDFnode::Sequence();
        for(const ASDFnode& item : items)
            list.Append(item);
        return list;
    }

    ASDFnode shiftNode(double offset)
    {
        ASDFnode node = transformNode("shift");
        node["offset"] = offset;
        return node;
    }

    ASDFnode scaleNode(double factor)
    {
        ASDFnode node = transformNode("scale");
        node["factor"] = factor;
        return node;
    }

    // A kind of transform the library doesn't know: factor * x^3
    class CubeTransform final : public ASDFtransform
    {
    public:
        explicit CubeTransform(double factor): ASDFtransform(1, 1), factor_(factor) {}
        void Evaluate(const double* in, double* out) const override { out[0] = factor_ * in[0] * in[0] * in[0]; }

    private:
        double factor_;
    };

    // Distance on the sky, degrees
    double skyDistance(double ra1, double dec1, double ra2, double dec2)
    {
        const double dra = std::fmod(ra1 - ra2 + 540., 360.) - 180.;
        return std::hypot(dra * std::cos(dec1 * std::numbers::pi / 180.), dec1 - dec2);
    }

    // The vertices (ra, dec) of an S_REGION "POLYGON ICRS ra1 dec1 ra2 dec2 ..."
    std::vector<std::pair<double,double>> polygon(const std::string& sRegion)
    {
        std::istringstream in(sRegion);
        in.imbue(std::locale::classic());
        std::string word;
        in >> word >> word;
        std::vector<std::pair<double,double>> vertices;
        double ra = 0., dec = 0.;
        while(in >> ra >> dec)
            vertices.emplace_back(ra, dec);
        return vertices;
    }

    // Largest distance between the vertices of an S_REGION and the nearest corner of the bounding box, for the
    // pipeline's way of taking the corners that fits best:
    //  - the edges of the box, as gwcs footprint();
    //  - the centres of the corner pixels;
    //  - the edges rounded to whole pixels, floor(edge + 0.5), as gwcs footprint(center=True) used by older
    //    romancal: (-0.5, 4087.5) gives (0, 4088), half a pixel off the two others on each axis.
    double sRegionError(const ASDFgwcs& gwcs, const std::string& sRegion)
    {
        const double inf = std::numeric_limits<double>::infinity();
        const auto vertices = polygon(sRegion);
        const auto& box = gwcs.BoundingBox();
        if(vertices.size() != 4 || box.size() != 2)
            return inf;

        using Edges = std::function<std::pair<double,double>(const std::pair<double,double>&)>;
        const Edges conventions[] = {
            [](const std::pair<double,double>& e) { return e; },
            [](const std::pair<double,double>& e) { return std::make_pair(e.first + 0.5, e.second - 0.5); },
            [](const std::pair<double,double>& e) { return std::make_pair(std::floor(e.first + 0.5), std::floor(e.second + 0.5)); }};

        double best = inf;
        for(const Edges& edges : conventions)
        {
            const auto [x0, x1] = edges(box[0]);
            const auto [y0, y1] = edges(box[1]);
            std::vector<std::vector<double>> corners;
            for(double x : {x0, x1})
                for(double y : {y0, y1})
                    corners.push_back(gwcs.PixelToWorld({x, y}, false));
            double worst = 0.;
            for(const auto& [ra, dec] : vertices)
            {
                double nearest = inf;
                for(const auto& corner : corners)
                    nearest = std::min(nearest, skyDistance(ra, dec, corner[0], corner[1]));
                worst = std::max(worst, nearest);
            }
            best = std::min(best, worst);
        }
        return best;
    }

    // Files downloaded by cmake/FetchSamples.cmake
    std::string sample(const std::string& name)
    {
        return std::string(ASDF2FITS_SAMPLES) + "/" + name;
    }
}

#define SKIP_WITHOUT(path) \
    do { if(!std::filesystem::exists(path)) GTEST_SKIP() << (path) << " is missing: cmake --build build --target fetch_samples"; } while(0)

TEST(ASDFtransform, KindIsTheTagWithoutItsVersion)
{
    EXPECT_EQ(ASDFtransform::KindOf("tag:stsci.edu:asdf/transform/shift-1.3.0"), "transform/shift");
    EXPECT_EQ(ASDFtransform::KindOf("asdf://asdf-format.org/transform/tags/shift-1.4.0"), "transform/shift");
    EXPECT_EQ(ASDFtransform::KindOf("tag:stsci.edu:gwcs/spherical_cartesian-1.2.0"), "gwcs/spherical_cartesian");
    EXPECT_EQ(ASDFtransform::KindOf(""), "");
    for(const char* kind : {"transform/compose", "transform/concatenate", "transform/shift", "transform/scale",
                            "transform/remap_axes", "transform/polynomial", "transform/affine", "transform/rotate_sequence_3d",
                            "transform/divide", "transform/constant", "transform/identity", "gwcs/spherical_cartesian",
                            "transform/gnomonic", "transform/rotate3d"})
        EXPECT_TRUE(ASDFtransform::IsKnown(kind)) << kind;
}

TEST(ASDFtransform, SimpleTransformsAndTheirInverses)
{
    // (x + 2.5) * 4, and back
    ASDFnode compose = transformNode("compose");
    compose["forward"] = nodeList({shiftNode(2.5), scaleNode(4.)});
    auto chain = ASDFtransform::Read(compose);
    EXPECT_EQ(chain->Kind(), "transform/compose");
    EXPECT_EQ(chain->Inputs(), 1u);
    EXPECT_DOUBLE_EQ((*chain)({1.})[0], 14.);
    ASSERT_TRUE(chain->Inverse());
    EXPECT_DOUBLE_EQ((*chain->Inverse())({14.})[0], 1.);

    // Mappings give each input back from the first output it went to
    ASDFnode remap = transformNode("remap_axes", "1.4.0");
    remap["mapping"] = numberList({1, 0, 1});
    auto mapping = ASDFtransform::Read(remap);
    EXPECT_EQ((*mapping)({7., 8.}), (std::vector<double>{8., 7., 8.}));
    ASSERT_TRUE(mapping->Inverse());
    EXPECT_EQ((*mapping->Inverse())({8., 7., 8.}), (std::vector<double>{7., 8.}));

    // Affine transforms are inverted as matrices
    ASDFnode affineNode = transformNode("affine", "1.4.0");
    affineNode["matrix"] = nodeList({numberList({1., 2.}), numberList({3., 4.})});
    affineNode["translation"] = numberList({5., 6.});
    auto affine = ASDFtransform::Read(affineNode);
    EXPECT_EQ((*affine)({1., 1.}), (std::vector<double>{8., 13.}));
    const std::vector<double> back = (*affine->Inverse())({8., 13.});
    EXPECT_NEAR(back[0], 1., 1e-12);
    EXPECT_NEAR(back[1], 1., 1e-12);

    // Rotations turn the frame, as astropy's: by +90 degrees about z, the x axis becomes -y
    ASDFnode quarter = transformNode("rotate_sequence_3d", "1.1.0");
    quarter["angles"] = numberList({90.});
    quarter["axes_order"] = ASDFnode("z");
    const std::vector<double> turned = (*ASDFtransform::Read(quarter))({1., 0., 0.});
    EXPECT_NEAR(turned[0], 0., 1e-15);
    EXPECT_NEAR(turned[1], -1., 1e-15);
    EXPECT_NEAR(turned[2], 0., 1e-15);

    ASDFnode rotation = transformNode("rotate_sequence_3d", "1.1.0");
    rotation["angles"] = numberList({10., 20., 30.});
    rotation["axes_order"] = ASDFnode("zyx");
    rotation["rotation_type"] = ASDFnode("cartesian");
    auto rotate = ASDFtransform::Read(rotation);
    const std::vector<double> rotated = (*rotate)({0.6, 0., 0.8});
    EXPECT_NEAR(std::hypot(rotated[0], rotated[1], rotated[2]), 1., 1e-15);
    const std::vector<double> unrotated = (*rotate->Inverse())(rotated);
    EXPECT_NEAR(unrotated[0], 0.6, 1e-15);
    EXPECT_NEAR(unrotated[1], 0., 1e-15);
    EXPECT_NEAR(unrotated[2], 0.8, 1e-15);

    // Longitude and latitude to a vector and back, the longitude in [0, 360)
    ASDFnode s2c = ASDFnode::Mapping("tag:stsci.edu:gwcs/spherical_cartesian-1.2.0");
    s2c["transform_type"] = ASDFnode("spherical_to_cartesian");
    s2c["wrap_lon_at"] = 360;
    auto toVector = ASDFtransform::Read(s2c);
    const std::vector<double> y = (*toVector)({90., 0.});
    EXPECT_NEAR(y[0], 0., 1e-15);
    EXPECT_NEAR(y[1], 1., 1e-15);
    EXPECT_NEAR(y[2], 0., 1e-15);
    const std::vector<double> lonLat = (*toVector->Inverse())({0., -1., 0.});
    EXPECT_NEAR(lonLat[0], 270., 1e-12);
    EXPECT_NEAR(lonLat[1], 0., 1e-12);
}

TEST(ASDFtransform, PolynomialsAndTheInverseGivenByTheFile)
{
    // The coefficient of x^i y^j is at [i, j]: 1 + 3x + 2y
    ASDFnode plane = transformNode("polynomial", "1.2.0");
    plane["coefficients"] = nodeList({numberList({1., 2., 0.}), numberList({3., 0., 0.}), numberList({0., 0., 0.})});
    auto p = ASDFtransform::Read(plane);
    EXPECT_EQ(p->Inputs(), 2u);
    EXPECT_DOUBLE_EQ((*p)({2., 5.})[0], 17.);
    EXPECT_FALSE(p->Inverse());

    // The domain is mapped to the window first: [0, 10] to [-1, 1]
    ASDFnode mapped = transformNode("polynomial", "1.2.0");
    mapped["coefficients"] = numberList({0., 1.});
    mapped["domain"] = numberList({0., 10.});
    mapped["window"] = numberList({-1., 1.});
    auto m = ASDFtransform::Read(mapped);
    // Mapped as numpy and astropy do, offset + scale*x: 0.2 isn't exact, so the middle comes out at about
    // 1e-17, not 0, as in astropy. A relative comparison (EXPECT_DOUBLE_EQ) can't be used next to 0.
    EXPECT_NEAR((*m)({5.})[0], 0., 1e-15);
    EXPECT_NEAR((*m)({0.})[0], -1., 1e-15);
    EXPECT_NEAR((*m)({10.})[0], 1., 1e-15);

    // 1 + 2x, with the inverse the file gives: (y - 1) / 2
    ASDFnode line = transformNode("polynomial", "1.2.0");
    line["coefficients"] = numberList({1., 2.});
    ASDFnode inverse = transformNode("compose");
    inverse["forward"] = nodeList({shiftNode(-1.), scaleNode(0.5)});
    line["inverse"] = inverse;
    auto l = ASDFtransform::Read(line);
    EXPECT_DOUBLE_EQ((*l)({2.})[0], 5.);
    ASSERT_TRUE(l->Inverse());
    EXPECT_DOUBLE_EQ((*l->Inverse())({5.})[0], 2.);
}

TEST(ASDFtransform, UnknownKindsFailUntilRegistered)
{
    ASDFnode node = ASDFnode::Mapping("tag:example.org:test/transform/cube-1.0.0");
    node["factor"] = 2.;
    ASSERT_FALSE(ASDFtransform::IsKnown("transform/cube"));
    EXPECT_THROW(ASDFtransform::Read(node), WCSexception);

    ASDFtransform::Register("transform/cube", [](const ASDFnode& n)
    {
        return std::make_shared<CubeTransform>(ASDFtransform::Number(n, "factor", 1.));
    });
    EXPECT_TRUE(ASDFtransform::IsKnown("transform/cube"));

    auto cube = ASDFtransform::Read(node);
    EXPECT_EQ(cube->Kind(), "transform/cube");
    EXPECT_DOUBLE_EQ((*cube)({3.})[0], 54.);
    EXPECT_FALSE(cube->Inverse());

    // Inside the transforms of the library too: 2 (x + 1)^3
    ASDFnode compose = transformNode("compose");
    compose["forward"] = nodeList({shiftNode(1.), node});
    EXPECT_DOUBLE_EQ((*ASDFtransform::Read(compose))({2.})[0], 54.);
    EXPECT_FALSE(ASDFtransform::Read(compose)->Inverse());
}

TEST(ASDFgwcs, RomanFixtureTakesItsReferencePointToTheSky)
{
    ASDFmanager file(fixture("roman_l2_wcs"));
    ASDFnode tree = file.GetTree();
    auto gwcs = ASDFgwcs::Read(tree.At("roman/meta/wcs"));

    EXPECT_EQ(gwcs->InputFrame().name, "detector");
    EXPECT_EQ(gwcs->OutputFrame().name, "world");
    EXPECT_TRUE(gwcs->OutputFrame().IsCelestial());
    EXPECT_EQ(gwcs->OutputFrame().referenceFrame, "icrs");
    EXPECT_EQ(gwcs->OutputFrame().units, (std::vector<std::string>{"deg", "deg"}));
    EXPECT_EQ(gwcs->PixelAxes(), 2u);
    EXPECT_EQ(gwcs->WorldAxes(), 2u);
    EXPECT_TRUE(gwcs->HasInverse());
    ASSERT_EQ(gwcs->BoundingBox().size(), 2u);
    EXPECT_EQ(gwcs->BoundingBox()[0], std::make_pair(-0.5, 4087.5));
    EXPECT_EQ(gwcs->BoundingBox()[1], std::make_pair(-0.5, 4087.5));

    // The pipeline points (v2_ref, v3_ref) of the telescope at (ra_ref, dec_ref)
    const ASDFnode& info = tree.At("roman/meta/wcsinfo");
    const double v2 = info.At("v2_ref").Get<double>(), v3 = info.At("v3_ref").Get<double>();
    const double ra = info.At("ra_ref").Get<double>(), dec = info.At("dec_ref").Get<double>();
    const std::vector<double> sky = gwcs->Convert("v2v3corr", "world", {v2, v3});
    EXPECT_LT(skyDistance(sky[0], sky[1], ra, dec), 1e-10) << sky[0] << " " << sky[1];
    const std::vector<double> v2v3 = gwcs->Convert("world", "v2v3corr", {ra, dec});
    EXPECT_NEAR(v2v3[0], v2, 1e-6);
    EXPECT_NEAR(v2v3[1], v3, 1e-6);
}

TEST(ASDFgwcs, RomanFixtureMatchesTheCornersOfItsSRegion)
{
    ASDFmanager file(fixture("roman_l2_wcs"));
    ASDFnode tree = file.GetTree();
    auto gwcs = ASDFgwcs::Read(tree.At("roman/meta/wcs"));

    // S_REGION is written by the pipeline with 9 decimals: 1e-8 degree is 36 micro-arcseconds
    const double error = sRegionError(*gwcs, tree.At("roman/meta/wcsinfo/s_region").Get<std::string>());
    EXPECT_LT(error, 1e-8) << "largest distance to a corner: " << error * 3600. << " arcsec";
}

TEST(ASDFgwcs, RomanFixtureRoundTripsInsideTheBoundingBoxOnly)
{
    ASDFmanager file(fixture("roman_l2_wcs"));
    auto gwcs = ASDFgwcs::Read(file.GetTree().At("roman/meta/wcs"));

    for(const std::vector<double>& pixel : std::vector<std::vector<double>>{{0., 0.}, {2043.5, 2043.5}, {4087., 0.}, {100.25, 3900.75}, {4087., 4087.}})
    {
        const std::vector<double> world = gwcs->PixelToWorld(pixel);
        ASSERT_TRUE(std::isfinite(world[0]) && std::isfinite(world[1]));
        const std::vector<double> back = gwcs->WorldToPixel(world);
        EXPECT_NEAR(back[0], pixel[0], 0.01) << pixel[0] << "," << pixel[1];
        EXPECT_NEAR(back[1], pixel[1], 0.01) << pixel[0] << "," << pixel[1];
    }

    // Outside the bounding box: NaN, unless asked otherwise
    EXPECT_TRUE(std::isnan(gwcs->PixelToWorld({-0.6, 10.})[0]));
    EXPECT_TRUE(std::isfinite(gwcs->PixelToWorld({-0.6, 10.}, false)[0]));
    const std::vector<double> far = gwcs->PixelToWorld({4500., 4500.}, false);
    EXPECT_TRUE(std::isnan(gwcs->WorldToPixel(far)[0]));
}

TEST(ASDFgwcs, PixelFoundByNewtonWhenAStepHasNoInverse)
{
    TempFile path(".asdf");
    writeRoman(path.path(), kPolynomialWcs);
    ASDFmanager file(path.path());
    auto gwcs = ASDFgwcs::Read(file.GetTree().At("roman/meta/wcs"));
    EXPECT_FALSE(gwcs->HasInverse());
    EXPECT_TRUE(gwcs->BoundingBox().empty());

    const std::vector<double> world = gwcs->PixelToWorld({10., 20.});
    EXPECT_NEAR(world[0], 6.03, 1e-12);
    EXPECT_NEAR(world[1], 10.06, 1e-12);

    for(const std::vector<double>& pixel : std::vector<std::vector<double>>{{10., 20.}, {-3.5, 61.25}, {40., 2.}})
    {
        const std::vector<double> back = gwcs->WorldToPixel(gwcs->PixelToWorld(pixel));
        EXPECT_NEAR(back[0], pixel[0], 1e-8);
        EXPECT_NEAR(back[1], pixel[1], 1e-8);
    }
}

TEST(ASDFimage, WorldCoordinatesComeFromTheGWCS)
{
    TempFile path(".asdf");
    writeRoman(path.path(), kLinearWcs);
    ASDFmanager file(path.path());
    std::shared_ptr<FITScube> image = file.GetImage("roman.data");
    auto asdf = asASDF(image);
    ASSERT_TRUE(asdf->GWCS());
    EXPECT_EQ(asdf->GWCSPath(), "roman/meta/wcs");
    EXPECT_EQ(asdf->GWCS()->OutputFrame().axesNames, (std::vector<std::string>{"u", "v"}));

    // What the axes are, through FITScube
    EXPECT_TRUE(image->HasWorldCoordinates());
    EXPECT_EQ(image->getNumberOfWCS(), 0u);     // no FITS WCS in the header
    EXPECT_EQ(image->WorldAxisNames(), (std::vector<std::string>{"u", "v"}));
    EXPECT_EQ(image->WorldAxisUnits(), (std::vector<std::string>{"deg", "deg"}));

    // Every form of WorldCoordinates() of FITScube comes to the GWCS
    const worldCoords w = image->WorldCoordinates(pixelCoords{10., 20.});
    ASSERT_EQ(w.size(), 2u);
    EXPECT_NEAR(w[0], -0.0215, 1e-15);
    EXPECT_NEAR(w[1], -0.007, 1e-15);
    EXPECT_EQ(image->WorldCoordinates(std::vector<size_t>{10, 20}), w);
    EXPECT_EQ(image->WorldCoordinates(size_t(20 * kNX + 10)), w);
    const worldVectors several = image->WorldCoordinatesVector(pixelVectors{{10., 20.}, {0., 0.}});
    ASSERT_EQ(several.size(), 2u);
    EXPECT_EQ(several[0], w);
    EXPECT_NEAR(several[1][0], -0.0315, 1e-15);

    const pixelCoords p = image->World2Pixel(w);
    EXPECT_NEAR(p[0], 10., 1e-9);
    EXPECT_NEAR(p[1], 20., 1e-9);
    const pixelVectors ps = image->World2PixelVector(worldVectors{w, several[1]});
    EXPECT_NEAR(ps[1][0], 0., 1e-9);
    EXPECT_NEAR(ps[1][1], 0., 1e-9);

    // Outside the bounding box
    EXPECT_TRUE(std::isnan(image->WorldCoordinates(pixelCoords{-1., 0.})[0]));
    EXPECT_TRUE(std::isnan(image->World2Pixel(worldCoords{1., 1.})[0]));
}

TEST(ASDFimage, AnUnusableGWCSLeavesTheFITSWCS)
{
    TempFile path(".asdf");
    writeRoman(path.path());
    ASDFmanager file(path.path());
    std::shared_ptr<FITScube> image = file.GetImage("roman.data");
    auto asdf = asASDF(image);
    EXPECT_FALSE(asdf->GWCS());
    EXPECT_EQ(asdf->GWCSPath(), "roman/meta/wcs");
    EXPECT_FALSE(image->HasWorldCoordinates());
    EXPECT_TRUE(image->WorldAxisNames().empty());
    // No WCS at all: FITScube gives the pixel back
    EXPECT_EQ(image->WorldCoordinates(pixelCoords{3., 4.}), (worldCoords{3., 4.}));
}

TEST(ASDFimage, WindowShiftsTheGWCS)
{
    TempFile in(".asdf"), out(".asdf");
    writeRoman(in.path(), kLinearWcs);
    ASDFmanager file(in.path());
    std::shared_ptr<FITScube> image = file.GetImage("roman.data");
    std::shared_ptr<FITScube> cut = image->Window(10, 20, 30, 16);
    auto asdf = asASDF(cut);
    ASSERT_TRUE(asdf->GWCS());

    for(const auto& [x, y] : std::vector<std::pair<double,double>>{{0., 0.}, {29., 15.}, {12.5, 3.25}})
    {
        const worldCoords w = cut->WorldCoordinates(pixelCoords{x, y});
        const worldCoords e = image->WorldCoordinates(pixelCoords{x + 10., y + 20.});
        EXPECT_NEAR(w[0], e[0], 1e-15) << x << "," << y;
        EXPECT_NEAR(w[1], e[1], 1e-15) << x << "," << y;
    }
    const pixelCoords p = cut->World2Pixel(image->WorldCoordinates(pixelCoords{12., 25.}));
    EXPECT_NEAR(p[0], 2., 1e-9);
    EXPECT_NEAR(p[1], 5., 1e-9);

    // The bounding box is the window's
    EXPECT_EQ(asdf->GWCS()->BoundingBox()[0], std::make_pair(-0.5, 29.5));
    EXPECT_EQ(asdf->GWCS()->BoundingBox()[1], std::make_pair(-0.5, 15.5));
    EXPECT_TRUE(std::isnan(cut->WorldCoordinates(pixelCoords{30., 0.})[0]));
    ASSERT_EQ(asdf->PendingHistory().size(), 1u);
    EXPECT_NE(asdf->PendingHistory()[0].find("GWCS shifted to the window: roman/meta/wcs"), std::string::npos);
    EXPECT_EQ(asdf->Meta().At("roman/meta/wcs/steps/0/transform").Tag(), "tag:stsci.edu:asdf/transform/compose-1.3.0");

    // The image it was cut from is unchanged
    EXPECT_EQ(asASDF(image)->GWCS()->BoundingBox()[0], std::make_pair(-0.5, 63.5));

    // The shifted GWCS is written, and reads back the same
    asdf->WriteASDF(out.path());
    ASDFmanager back(out.path());
    std::shared_ptr<FITScube> reread = back.GetImage("roman.data");
    const worldCoords w = reread->WorldCoordinates(pixelCoords{12.5, 3.25});
    const worldCoords e = image->WorldCoordinates(pixelCoords{22.5, 23.25});
    EXPECT_NEAR(w[0], e[0], 1e-15);
    EXPECT_NEAR(w[1], e[1], 1e-15);

    // A window read from the file gets the same GWCS
    std::shared_ptr<FITScube> read = file.GetImageWindow("roman.data", {10, 20}, {30, 16});
    EXPECT_EQ(read->WorldCoordinates(pixelCoords{12.5, 3.25}), cut->WorldCoordinates(pixelCoords{12.5, 3.25}));
}

// ---------------------------------------------------------------------- Real files, downloaded on request

TEST(Samples, JwstNircamGWCS)
{
    const std::string path = sample("gwcs_nircamwcs.asdf");
    SKIP_WITHOUT(path);
    ASDFmanager file(path);
    auto gwcs = ASDFgwcs::Read(file.GetTree().At("wcs"));
    EXPECT_EQ(gwcs->InputFrame().name, "detector");
    EXPECT_EQ(gwcs->OutputFrame().name, "world");
    EXPECT_TRUE(gwcs->HasInverse());

    // The angles of its rotation to the sky: (v2_ref, -v3_ref, roll, dec_ref, -ra_ref)
    const double v2 = 0.023917627222222224 * 3600., v3 = -0.13700764222222223 * 3600.;
    const double ra = 5.868934166666667, dec = -71.99550858333333;
    const std::vector<double> sky = gwcs->Convert("v2v3", "world", {v2, v3});
    EXPECT_LT(skyDistance(sky[0], sky[1], ra, dec), 1e-10) << sky[0] << " " << sky[1];

    // The way back uses the inverse distortion polynomials of the file, fitted to the forward ones: ~0.02 pixel
    // at the edges, 0 at the centre (gwcs' own test_iter_inv removes them to test an exact round trip)
    for(const std::vector<double>& pixel : std::vector<std::vector<double>>{{0., 0.}, {1023.5, 1023.5}, {2047., 2047.}, {100., 1900.}})
    {
        const std::vector<double> back = gwcs->WorldToPixel(gwcs->PixelToWorld(pixel));
        EXPECT_NEAR(back[0], pixel[0], 0.05) << pixel[0] << "," << pixel[1];
        EXPECT_NEAR(back[1], pixel[1], 0.05) << pixel[0] << "," << pixel[1];
    }
}

TEST(Samples, JwstMiriGWCS)
{
    const std::string path = sample("gwcs_miriwcs.asdf");
    SKIP_WITHOUT(path);
    ASDFmanager file(path);
    // Some inverses inside the file don't fit their transform (identity(2) for the mapping (0, 1, 0, 1)): they
    // are ignored, the whole step has its own
    auto gwcs = ASDFgwcs::Read(file.GetTree().At("wcs"));
    EXPECT_TRUE(gwcs->OutputFrame().IsCelestial());
    EXPECT_TRUE(gwcs->HasInverse());
    ASSERT_EQ(gwcs->BoundingBox().size(), 2u);

    // As for NIRCam, the way back uses fitted inverse polynomials
    const auto& box = gwcs->BoundingBox();
    for(double fx : {0.1, 0.5, 0.9})
        for(double fy : {0.1, 0.5, 0.9})
        {
            const std::vector<double> pixel{box[0].first + fx * (box[0].second - box[0].first),
                                            box[1].first + fy * (box[1].second - box[1].first)};
            const std::vector<double> world = gwcs->PixelToWorld(pixel);
            ASSERT_TRUE(std::isfinite(world[0]) && std::isfinite(world[1]));
            const std::vector<double> back = gwcs->WorldToPixel(world);
            EXPECT_NEAR(back[0], pixel[0], 0.05) << pixel[0] << "," << pixel[1];
            EXPECT_NEAR(back[1], pixel[1], 0.05) << pixel[0] << "," << pixel[1];
        }
}

TEST(Samples, RomanL2ImageWindowKeepsItsSky)
{
    const std::string path = sample("roman_wfi01_f106_cal.asdf");
    SKIP_WITHOUT(path);
    ASDFmanager file(path);
    ASDFnode tree = file.GetTree();
    auto full = ASDFgwcs::Read(tree.At("roman/meta/wcs"));
    const double error = sRegionError(*full, tree.At("roman/meta/wcsinfo/s_region").Get<std::string>());
    EXPECT_LT(error, 1e-8) << "largest distance to a corner: " << error * 3600. << " arcsec";

    std::shared_ptr<FITScube> cut = file.GetImageWindow("roman.data", {1000, 2000}, {64, 32});
    ASSERT_TRUE(asASDF(cut)->GWCS());
    for(const auto& [x, y] : std::vector<std::pair<double,double>>{{0., 0.}, {63., 31.}, {10.5, 20.25}})
    {
        const worldCoords w = cut->WorldCoordinates(pixelCoords{x, y});
        const std::vector<double> e = full->PixelToWorld({x + 1000., y + 2000.});
        EXPECT_NEAR(w[0], e[0], 1e-12);
        EXPECT_NEAR(w[1], e[1], 1e-12);
        const pixelCoords p = cut->World2Pixel(w);
        EXPECT_NEAR(p[0], x, 0.01);
        EXPECT_NEAR(p[1], y, 0.01);
    }
}

// ---------------------------------------------------------------------- Known bugs of libasdf, not patched
//
// Each test checks that a bug of libasdf is still there. When one fails, the libasdf now built has the bug
// fixed: remove the workaround it names, the test itself, and its line in LIBASDF_KNOWN_BUGS (CMakeLists.txt).

// asdf_set_*(file, "a/b") fails when the mapping "a" is already in the tree: asdf_node_materialize_path
// (libasdf/src/value.c) appends the existing "a" to its parent again, and reports that as out of memory,
// without an error message. The key is lost. Workaround: writeMultiHDU() builds each mapping whole.
TEST(LibasdfKnownBugs, SetUnderAnExistingMapping)
{
    asdf_file_t* file = asdf_open_mem(nullptr, 0);
    ASSERT_NE(file, nullptr);
    EXPECT_EQ(asdf_set_string0(file, "a/x", "one"), ASDF_VALUE_OK);     // creates a
    const asdf_value_err_t second = asdf_set_string0(file, "a/y", "two");
    asdf_close(file);

    EXPECT_EQ(second, ASDF_VALUE_ERR_OOM)
        << "libasdf " LIBASDF_VERSION " now adds a key to an existing mapping (returned " << second << "): the bug is "
           "fixed. writeMultiHDU() may use asdf_set_*(file, \"group/key\") again; remove this test, the comment above "
           "setNdarray() and the entry of LIBASDF_KNOWN_BUGS in CMakeLists.txt.";
}
