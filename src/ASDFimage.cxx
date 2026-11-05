//
//  ASDFimage.cxx
//  asdf2fits
//
//  Licensed under PolyForm Noncommercial 1.0.0
//  You may use, share and adapt this code for noncommercial purposes.
//  Commercial use requires a separate licence from the author.
//  Licence text: https://polyformproject.org/licenses/noncommercial/1.0.0/

#include "ASDFimage.h"

#include <bit>
#include <iterator>

#include <cmath>
#include <limits>
#include <type_traits>

#include "ASDFheader.h"

namespace DSL::asdfdetail
{
    namespace
    {
        const char* const kWho = "asdf2fits";

        template<typename T> asdf_scalar_datatype_t scalarType()
        {
            if constexpr (std::is_same_v<T, int8_t>)   return ASDF_DATATYPE_INT8;
            if constexpr (std::is_same_v<T, uint8_t>)  return ASDF_DATATYPE_UINT8;
            if constexpr (std::is_same_v<T, int16_t>)  return ASDF_DATATYPE_INT16;
            if constexpr (std::is_same_v<T, uint16_t>) return ASDF_DATATYPE_UINT16;
            if constexpr (std::is_same_v<T, int32_t>)  return ASDF_DATATYPE_INT32;
            if constexpr (std::is_same_v<T, uint32_t>) return ASDF_DATATYPE_UINT32;
            if constexpr (std::is_integral_v<T> && std::is_signed_v<T>)   return ASDF_DATATYPE_INT64;
            if constexpr (std::is_integral_v<T> && std::is_unsigned_v<T>) return ASDF_DATATYPE_UINT64;
            if constexpr (std::is_same_v<T, float>)    return ASDF_DATATYPE_FLOAT32;
            if constexpr (std::is_same_v<T, double>)   return ASDF_DATATYPE_FLOAT64;
        }

        template<typename T>
        std::shared_ptr<FITScube> load(asdf_ndarray_t* ndarray, const FITShdu& header)
        {
            const std::vector<uint64_t> axes = fitsAxes(ndarray);
            auto image = std::make_shared<FITSimg<T>>(std::vector<size_t>(axes.begin(), axes.end()));

            void* dst = image->template raw_data<T>();
            const asdf_ndarray_err_t err = asdf_ndarray_read_all(ndarray, scalarType<T>(), &dst);
            if(err != ASDF_NDARRAY_OK)
                throw FITSexception(READ_ERROR, kWho, "readImage", "libasdf error "+std::to_string(err)+" reading the array");

            image->HDU().Import(header, {"BLANK"});

            for(const char* key : {"COMMENT", "HISTORY"})
                if(header.Exists(key))
                    image->HDU().ValueForKey(key, header.GetEntry(key)->second.value(), fChar);
            if(header.Exists("EXTNAME"))
                image->SetName(header.GetValueForKey("EXTNAME"));

            // Masked pixels are recognised as FITSimg does when it reads a FITS file
            const std::valarray<T>& pixels = *image->template GetData<T>();
            std::valarray<bool> masked(false, pixels.size());
            bool any = false;
            if constexpr (std::is_floating_point_v<T>)
            {
                for(size_t i = 0; i < pixels.size(); i++)
                    if(std::isnan(pixels[i]))
                        any = masked[i] = true;
            }
            else if(header.Exists("BLANK"))
            {
                const T blank = std::is_unsigned_v<T> ? static_cast<T>(header.GetUInt64ValueForKey("BLANK"))
                                                      : static_cast<T>(header.GetInt64ValueForKey("BLANK"));
                image->Blank(blank);
                for(size_t i = 0; i < pixels.size(); i++)
                    if(pixels[i] == blank)
                        any = masked[i] = true;
            }
            if(any)
                image->MaskPixels(masked);

            if(header.Exists("CTYPE1") || header.Exists("CRPIX1") || header.Exists("WCSAXES"))
                image->reLoadWCS();

            return image;
        }

        template<typename T>
        bool store(asdf_file_t* file, const FITScube& cube, FITShdu& header, asdf_value_t*& out)
        {
            const std::valarray<T>* pixels = cube.GetData<T>();
            if(!pixels)
                return false;

            std::vector<uint64_t> shape;
            std::vector<uint64_t> naxis;
            for(size_t i = cube.GetDimension(); i > 0; i--)
                shape.push_back(cube.Size(i));
            for(size_t i = 1; i <= cube.GetDimension(); i++)
                naxis.push_back(cube.Size(i));

            asdf_ndarray_t ndarray{};
            ndarray.ndim = static_cast<uint32_t>(shape.size());
            ndarray.shape = shape.data();
            ndarray.datatype.type = scalarType<T>();
            ndarray.byteorder = (std::endian::native == std::endian::little) ? ASDF_BYTEORDER_LITTLE : ASDF_BYTEORDER_BIG;

            T* buffer = static_cast<T*>(asdf_ndarray_data_alloc(&ndarray));
            if(!buffer)
                throw FITSexception(MEMORY_ALLOCATION, kWho, "writeImage", lastError(file));

            setImageStructure(header, scalarType<T>(), naxis);
            if(header.Exists("BLANK"))
                header.DeleteKey("BLANK");

            const std::valarray<bool> mask = cube.GetMask();
            const bool anyMasked = mask.size() == pixels->size() && mask.max();

            T blank = T();
            if constexpr (std::is_floating_point_v<T>)
                blank = std::numeric_limits<T>::quiet_NaN();
            else
            {
                blank = std::is_unsigned_v<T> ? std::numeric_limits<T>::max() : std::numeric_limits<T>::min();
                if(anyMasked)
                {
                    if constexpr (std::is_unsigned_v<T>)
                        header.ValueForKey("BLANK", static_cast<uint64_t>(blank));
                    else
                        header.ValueForKey("BLANK", static_cast<int64_t>(blank));
                }
            }

            for(size_t i = 0; i < pixels->size(); i++)
                buffer[i] = (anyMasked && mask[i]) ? blank : (*pixels)[i];

            out = asdf_value_of_ndarray(file, &ndarray);
            if(!out)
            {
                asdf_ndarray_data_dealloc(&ndarray);
                throw FITSexception(WRITE_ERROR, kWho, "writeImage", lastError(file));
            }
            return true;
        }
    }

    std::vector<uint64_t> fitsAxes(const asdf_ndarray_t* ndarray)
    {
        return std::vector<uint64_t>(std::make_reverse_iterator(ndarray->shape + ndarray->ndim),
                                     std::make_reverse_iterator(ndarray->shape));
    }

    bool isImageType(asdf_scalar_datatype_t type)
    {
        switch(type)
        {
            case ASDF_DATATYPE_INT8:  case ASDF_DATATYPE_UINT8:
            case ASDF_DATATYPE_INT16: case ASDF_DATATYPE_UINT16:
            case ASDF_DATATYPE_INT32: case ASDF_DATATYPE_UINT32:
            case ASDF_DATATYPE_INT64: case ASDF_DATATYPE_UINT64:
            case ASDF_DATATYPE_FLOAT16: case ASDF_DATATYPE_FLOAT32: case ASDF_DATATYPE_FLOAT64:
            case ASDF_DATATYPE_BOOL8:
                return true;
            default:
                return false;
        }
    }

    std::shared_ptr<FITScube> readImage(asdf_ndarray_t* ndarray, const FITShdu& header)
    {
        switch(ndarray->datatype.type)
        {
            case ASDF_DATATYPE_INT8:    return load<int8_t>(ndarray, header);
            case ASDF_DATATYPE_BOOL8:
            case ASDF_DATATYPE_UINT8:   return load<uint8_t>(ndarray, header);
            case ASDF_DATATYPE_INT16:   return load<int16_t>(ndarray, header);
            case ASDF_DATATYPE_UINT16:  return load<uint16_t>(ndarray, header);
            case ASDF_DATATYPE_INT32:   return load<int32_t>(ndarray, header);
            case ASDF_DATATYPE_UINT32:  return load<uint32_t>(ndarray, header);
            case ASDF_DATATYPE_INT64:   return load<int64_t>(ndarray, header);
            case ASDF_DATATYPE_UINT64:  return load<uint64_t>(ndarray, header);
            case ASDF_DATATYPE_FLOAT16:
            case ASDF_DATATYPE_FLOAT32: return load<float>(ndarray, header);
            case ASDF_DATATYPE_FLOAT64: return load<double>(ndarray, header);
            default:
                throw FITSexception(BAD_DATATYPE, kWho, "readImage", std::string("no FITS image holds datatype ")+asdf_scalar_datatype_to_string(ndarray->datatype.type));
        }
    }

    asdf_value_t* writeImage(asdf_file_t* file, const FITScube& cube, FITShdu& header)
    {
        asdf_value_t* out = nullptr;
        const bool done = store<int8_t>(file, cube, header, out)   || store<uint8_t>(file, cube, header, out)
                       || store<int16_t>(file, cube, header, out)  || store<uint16_t>(file, cube, header, out)
                       || store<int32_t>(file, cube, header, out)  || store<uint32_t>(file, cube, header, out)
                       || store<int64_t>(file, cube, header, out)  || store<uint64_t>(file, cube, header, out)
#if defined(__APPLE__)
                       || store<size_t>(file, cube, header, out)
#endif
                       || store<float>(file, cube, header, out)    || store<double>(file, cube, header, out);
        if(!done)
            throw FITSexception(BAD_DATATYPE, kWho, "writeImage", "the cube holds no data of a known type");
        return out;
    }
}
