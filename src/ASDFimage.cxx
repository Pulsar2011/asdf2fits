//
//  ASDFimage.cxx
//  asdf2fits
//
//  Written by Claude (Anthropic's AI model) for William Gillard.
//  Copyright (c) 2026 William Gillard.
//
//  SPDX-License-Identifier: BSD-3-Clause
//  Licence text: see the LICENSE file at the root of this project.

/**
 * @file ASDFimage.cxx
 * @brief Implementation of the conversion between ndarrays and FITSimg\<T\>.
 */
#include "ASDFimage.h"

#include <bit>
#include <iterator>

#include <cmath>
#include <limits>
#include <regex>
#include <type_traits>

#include "ASDFcube.h"
#include "ASDFheader.h"

namespace DSL::asdfdetail
{
    namespace
    {
        const char* const kWho = "asdf2fits";   //!< Name given in the exceptions thrown here

        /**
         * @brief Get the ndarray element type of a pixel type.
         * @tparam T Pixel type: any integer type, float or double.
         * @return The matching libasdf type.
         */
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

        /**
         * @brief Get the value given to masked pixels.
         * @tparam T Pixel type.
         * @return NaN for a real type; the lowest value of a signed type, the highest of an unsigned one.
         */
        template<typename T> T blankValue()
        {
            if constexpr (std::is_floating_point_v<T>)
                return std::numeric_limits<T>::quiet_NaN();
            else
                return std::is_unsigned_v<T> ? std::numeric_limits<T>::max() : std::numeric_limits<T>::min();
        }

        /**
         * @brief Check that a tile lies inside an ndarray.
         * @param ndarray Ndarray to read.
         * @param tile    Tile to read.
         * @throws FITSexception BAD_DIMEN if it doesn't.
         */
        void checkTile(const asdf_ndarray_t* ndarray, const Tile& tile)
        {
            if(tile.whole())
                return;
            if(tile.origin.size() != ndarray->ndim || tile.size.size() != ndarray->ndim || tile.dropped >= ndarray->ndim + 1)
                throw FITSexception(BAD_DIMEN, kWho, "readImage", "the tile doesn't have the "+std::to_string(ndarray->ndim)+" axes of the array");
            for(uint32_t i = 0; i < ndarray->ndim; i++)
                if(tile.origin[i] + tile.size[i] > ndarray->shape[i])
                    throw FITSexception(BAD_DIMEN, kWho, "readImage", "the tile goes past the end of axis "+std::to_string(ndarray->ndim - i)+" of the array");
        }

        /**
         * @brief Load an ndarray, or part of it, into a new image of pixel type @p T.
         * @details Masks the pixels that are NaN, or equal to BLANK in an integer image, and loads the WCS if
         * the header describes one.
         * @tparam T Pixel type of the image; the data are converted to it.
         * @param ndarray   Ndarray giving the pixels and the shape.
         * @param header    Header giving the keywords of the image.
         * @param tile      Part to read.
         * @param asdfImage Make an ASDFimg\<T\>.
         * @return The new image.
         * @throws FITSexception BAD_DIMEN if the tile isn't inside the array, READ_ERROR if reading fails.
         */
        template<typename T>
        std::shared_ptr<FITScube> load(asdf_ndarray_t* ndarray, const FITShdu& header, const Tile& tile, bool asdfImage)
        {
            checkTile(ndarray, tile);
            const std::vector<uint64_t> shape = tile.whole() ? std::vector<uint64_t>(ndarray->shape, ndarray->shape + ndarray->ndim) : tile.size;

            // FITS order, without the leading axes the tile leaves out
            std::vector<size_t> axes(shape.rbegin(), shape.rend() - static_cast<std::ptrdiff_t>(std::min(tile.dropped, shape.size())));
            if(axes.empty())
                axes.push_back(1);

            std::shared_ptr<FITSimg<T>> image;
            if(asdfImage)
                image = std::make_shared<ASDFimg<T>>(axes);
            else
                image = std::make_shared<FITSimg<T>>(axes);

            void* dst = image->template raw_data<T>();
            if(dst)
            {
                const asdf_ndarray_err_t err = tile.whole()
                    ? asdf_ndarray_read_all(ndarray, scalarType<T>(), &dst)
                    : asdf_ndarray_read_tile_ndim(ndarray, tile.origin.data(), tile.size.data(), scalarType<T>(), &dst);
                if(err != ASDF_NDARRAY_OK)
                    throw FITSexception(READ_ERROR, kWho, "readImage", "libasdf error "+std::to_string(err)+" reading the array");
            }

            if(tile.whole())
                image->HDU().Import(header, {"BLANK"});
            else
            {
                // The header describes the whole array: shift the reference pixel of its WCS to the tile
                FITShdu shifted(header);
                static const std::regex crpix("CRPIX([0-9]+)[A-Z]?");
                std::vector<std::string> keys;
                for(auto it = shifted.begin(); it != shifted.end(); ++it)
                    if(std::regex_match(it->first, crpix))
                        keys.push_back(it->first);
                for(const std::string& key : keys)
                {
                    const size_t axis = std::stoul(key.substr(5));
                    if(axis >= 1 && axis <= ndarray->ndim)
                        shifted.ValueForKey(key, shifted.GetDoubleValueForKey(key) - static_cast<double>(tile.origin[ndarray->ndim - axis]));
                }
                image->HDU().Import(shifted, {"BLANK"});

                // ... and give it the axes of the tile
                static const std::regex naxis("NAXIS[0-9]+");
                std::vector<std::string> stale;
                for(auto it = image->HDU().begin(); it != image->HDU().end(); ++it)
                    if(std::regex_match(it->first, naxis))
                        stale.push_back(it->first);
                for(const std::string& key : stale)
                    image->HDU().DeleteKey(key);
                image->HDU().ValueForKey("NAXIS", static_cast<int>(axes.size()));
                for(size_t i = 0; i < axes.size(); i++)
                    image->HDU().ValueForKey("NAXIS"+std::to_string(i + 1), axes[i]);
            }

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

        /**
         * @brief Copy a cube into a new ndarray value if its pixels are of type @p T.
         * @details Masked pixels take NaN or, in integer images, the value stored as BLANK.
         * @tparam T Pixel type to try.
         * @param file   File the value is created in.
         * @param cube   Cube to copy.
         * @param header Header of the output HDU; receives the structure keywords and BLANK if needed.
         * @param out    Receives the new value.
         * @return False if the cube holds another type, and nothing was done.
         * @throws FITSexception MEMORY_ALLOCATION if the buffer can't be allocated, WRITE_ERROR if the value can't be created.
         */
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
            const bool anyMasked = mask.size() == pixels->size() && mask.size() > 0 && mask.max();

            const T blank = blankValue<T>();
            if constexpr (!std::is_floating_point_v<T>)
            {
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

        /**
         * @brief Get the pixels of a cube if they are of type @p T.
         * @tparam T Pixel type to try.
         * @param cube Cube to write.
         * @param out  Receives the pixels.
         * @return False if the cube holds another type, and nothing was done.
         */
        template<typename T>
        bool pixelsOf(const FITScube& cube, PixelBlock& out)
        {
            const std::valarray<T>* pixels = cube.GetData<T>();
            if(!pixels)
                return false;

            out.type = scalarType<T>();
            out.shape.clear();
            for(size_t i = cube.GetDimension(); i > 0; i--)
                out.shape.push_back(cube.Size(i));
            out.bytes = pixels->size() * sizeof(T);

            const std::valarray<bool> mask = cube.GetMask();
            const bool anyMasked = mask.size() == pixels->size() && mask.size() > 0 && mask.max();
            if(!anyMasked)
            {
                out.data = pixels->size() ? &(*pixels)[0] : nullptr;
                return true;
            }

            const T blank = blankValue<T>();
            if constexpr (!std::is_floating_point_v<T>)
                out.blank = ASDFnode(blank);

            out.owned.resize(out.bytes);
            T* buffer = reinterpret_cast<T*>(out.owned.data());
            for(size_t i = 0; i < pixels->size(); i++)
                buffer[i] = mask[i] ? blank : (*pixels)[i];
            out.data = out.owned.data();
            return true;
        }

        /**
         * @brief Tell whether a cube holds pixels of type @p T.
         * @tparam T Pixel type to try.
         * @param cube The cube.
         * @param type Receives the libasdf type if it does.
         * @return True if it does.
         */
        template<typename T>
        bool typeOf(const FITScube& cube, asdf_scalar_datatype_t& type)
        {
            if(!cube.GetData<T>())
                return false;
            type = scalarType<T>();
            return true;
        }

        /**
         * @brief Copy a cube if it is a FITSimg\<T\>.
         * @tparam T Pixel type to try.
         * @param cube The cube.
         * @param out  Receives the copy.
         * @return False if the cube holds another type.
         */
        template<typename T>
        bool copyOf(const FITScube& cube, std::shared_ptr<FITScube>& out)
        {
            const FITSimg<T>* image = dynamic_cast<const FITSimg<T>*>(&cube);
            if(!image)
                return false;
            out = std::make_shared<FITSimg<T>>(*image);
            return true;
        }

        /**
         * @brief Copy part of a cube if it holds pixels of type @p T.
         * @tparam T Pixel type to try.
         * @param cube    The cube.
         * @param origin  First pixel, NAXIS1 first.
         * @param size    Length along each axis, NAXIS1 first.
         * @param dropped Last axes left out of the copy.
         * @param out     Receives the copy.
         * @return False if the cube holds another type.
         */
        template<typename T>
        bool windowOf(const FITScube& cube, const std::vector<size_t>& origin, const std::vector<size_t>& size, size_t dropped,
                      std::shared_ptr<FITScube>& out)
        {
            const std::valarray<T>* src = cube.GetData<T>();
            if(!src)
                return false;

            const size_t rank = cube.GetDimension();
            std::vector<size_t> axes(size.begin(), size.end() - static_cast<std::ptrdiff_t>(std::min(dropped, size.size())));
            if(axes.empty())
                axes.push_back(1);

            auto image = std::make_shared<FITSimg<T>>(axes);
            std::valarray<T>& dst = *image->template GetData<T>();
            const std::valarray<bool> srcMask = cube.GetMask();
            std::valarray<bool> mask(false, dst.size());
            bool any = false;

            // Walk the part, NAXIS1 fastest, and find each pixel in the cube
            std::vector<size_t> pos(rank, 0);
            for(size_t n = 0; n < dst.size(); n++)
            {
                size_t index = 0;
                size_t stride = 1;
                for(size_t k = 0; k < rank; k++)
                {
                    index += (origin[k] + pos[k]) * stride;
                    stride *= cube.Size(k + 1);
                }
                dst[n] = (*src)[index];
                if(index < srcMask.size() && srcMask[index])
                    any = mask[n] = true;

                for(size_t k = 0; k < rank; k++)
                {
                    if(++pos[k] < size[k])
                        break;
                    pos[k] = 0;
                }
            }

            if(any)
                image->MaskPixels(mask);
            if(cube.HDU().Exists("EXTNAME"))
                image->SetName(cube.HDU().GetValueForKey("EXTNAME"));
            out = image;
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

    std::shared_ptr<FITScube> readImage(asdf_ndarray_t* ndarray, const FITShdu& header, const Tile& tile, bool asdfImage)
    {
        switch(ndarray->datatype.type)
        {
            case ASDF_DATATYPE_INT8:    return load<int8_t>(ndarray, header, tile, asdfImage);
            case ASDF_DATATYPE_BOOL8:
            case ASDF_DATATYPE_UINT8:   return load<uint8_t>(ndarray, header, tile, asdfImage);
            case ASDF_DATATYPE_INT16:   return load<int16_t>(ndarray, header, tile, asdfImage);
            case ASDF_DATATYPE_UINT16:  return load<uint16_t>(ndarray, header, tile, asdfImage);
            case ASDF_DATATYPE_INT32:   return load<int32_t>(ndarray, header, tile, asdfImage);
            case ASDF_DATATYPE_UINT32:  return load<uint32_t>(ndarray, header, tile, asdfImage);
            case ASDF_DATATYPE_INT64:   return load<int64_t>(ndarray, header, tile, asdfImage);
            case ASDF_DATATYPE_UINT64:  return load<uint64_t>(ndarray, header, tile, asdfImage);
            case ASDF_DATATYPE_FLOAT16:
            case ASDF_DATATYPE_FLOAT32: return load<float>(ndarray, header, tile, asdfImage);
            case ASDF_DATATYPE_FLOAT64: return load<double>(ndarray, header, tile, asdfImage);
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

    PixelBlock cubePixels(const FITScube& cube)
    {
        PixelBlock out;
        const bool done = pixelsOf<int8_t>(cube, out)   || pixelsOf<uint8_t>(cube, out)
                       || pixelsOf<int16_t>(cube, out)  || pixelsOf<uint16_t>(cube, out)
                       || pixelsOf<int32_t>(cube, out)  || pixelsOf<uint32_t>(cube, out)
                       || pixelsOf<int64_t>(cube, out)  || pixelsOf<uint64_t>(cube, out)
#if defined(__APPLE__)
                       || pixelsOf<size_t>(cube, out)
#endif
                       || pixelsOf<float>(cube, out)    || pixelsOf<double>(cube, out);
        if(!done)
            throw FITSexception(BAD_DATATYPE, kWho, "cubePixels", "the cube holds no data of a known type");
        return out;
    }

    std::string datatypeName(asdf_scalar_datatype_t type)
    {
        const char* name = asdf_scalar_datatype_to_string(type);
        return name ? std::string(name) : std::string("unknown");
    }

    asdf_scalar_datatype_t cubeType(const FITScube& cube)
    {
        asdf_scalar_datatype_t type = ASDF_DATATYPE_UNKNOWN;
        [[maybe_unused]] const bool found =
            typeOf<int8_t>(cube, type)   || typeOf<uint8_t>(cube, type)  || typeOf<int16_t>(cube, type) || typeOf<uint16_t>(cube, type)
         || typeOf<int32_t>(cube, type)  || typeOf<uint32_t>(cube, type) || typeOf<int64_t>(cube, type) || typeOf<uint64_t>(cube, type)
#if defined(__APPLE__)
         || typeOf<size_t>(cube, type)
#endif
         || typeOf<float>(cube, type)    || typeOf<double>(cube, type);
        return type;
    }

    std::shared_ptr<FITScube> copyCube(const FITScube& cube)
    {
        std::shared_ptr<FITScube> out;
        const bool done = copyOf<int8_t>(cube, out)   || copyOf<uint8_t>(cube, out)  || copyOf<int16_t>(cube, out) || copyOf<uint16_t>(cube, out)
                       || copyOf<int32_t>(cube, out)  || copyOf<uint32_t>(cube, out) || copyOf<int64_t>(cube, out) || copyOf<uint64_t>(cube, out)
#if defined(__APPLE__)
                       || copyOf<size_t>(cube, out)
#endif
                       || copyOf<float>(cube, out)    || copyOf<double>(cube, out);
        if(!done)
            throw FITSexception(BAD_DATATYPE, kWho, "copyCube", "the cube holds no data of a known type");
        return out;
    }

    std::shared_ptr<FITScube> cubeWindow(const FITScube& cube, const std::vector<size_t>& origin, const std::vector<size_t>& size, size_t dropped)
    {
        std::shared_ptr<FITScube> out;
        const bool done = windowOf<int8_t>(cube, origin, size, dropped, out)   || windowOf<uint8_t>(cube, origin, size, dropped, out)
                       || windowOf<int16_t>(cube, origin, size, dropped, out)  || windowOf<uint16_t>(cube, origin, size, dropped, out)
                       || windowOf<int32_t>(cube, origin, size, dropped, out)  || windowOf<uint32_t>(cube, origin, size, dropped, out)
                       || windowOf<int64_t>(cube, origin, size, dropped, out)  || windowOf<uint64_t>(cube, origin, size, dropped, out)
#if defined(__APPLE__)
                       || windowOf<size_t>(cube, origin, size, dropped, out)
#endif
                       || windowOf<float>(cube, origin, size, dropped, out)    || windowOf<double>(cube, origin, size, dropped, out);
        if(!done)
            throw FITSexception(BAD_DATATYPE, kWho, "cubeWindow", "the cube holds no data of a known type");
        return out;
    }
}
