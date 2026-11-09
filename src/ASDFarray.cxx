//
//  ASDFarray.cxx
//  asdf2fits
//
//  Written by Claude (Anthropic's AI model) for William Gillard.
//  Copyright (c) 2026 William Gillard.
//
//  SPDX-License-Identifier: BSD-3-Clause
//  Licence text: see the LICENSE file at the root of this project.

/**
 * @file ASDFarray.cxx
 * @brief Implementation of DSL::ASDFarray and of the description of ndarrays in their file.
 */
#include "ASDFarray.h"

#include <charconv>
#include <numeric>

#include <DSTfits/FITSimg.h>

#include "ASDFimage.h"
#include "ASDFsource.h"

namespace DSL
{
    using namespace asdfdetail;

    namespace
    {
        const char* const kWho = "ASDFarray";   //!< Name given in the exceptions thrown here

        /**
         * @brief Bring a part given in the FITS order to the axes of a view.
         * @details Missing positions start at 0; missing or 0 lengths, and lengths past the end, go to the end.
         * @param axes   Axis lengths of the view, NAXIS1 first.
         * @param origin First pixel; receives one position per axis.
         * @param size   Lengths; receive one length per axis.
         * @param what   Name of the calling method, for the message.
         * @throws FITSexception BAD_DIMEN if there are more positions or lengths than axes, or the part starts outside.
         */
        void fitPart(const std::vector<size_t>& axes, std::vector<size_t>& origin, std::vector<size_t>& size, const char* what)
        {
            if(origin.size() > axes.size() || size.size() > axes.size())
                throw FITSexception(BAD_DIMEN, kWho, what, "the array has "+std::to_string(axes.size())+" axes");
            origin.resize(axes.size(), 0);
            size.resize(axes.size(), 0);
            for(size_t k = 0; k < axes.size(); k++)
            {
                if(origin[k] >= axes[k])
                    throw FITSexception(BAD_DIMEN, kWho, what, "position "+std::to_string(origin[k])+" is past the end of axis "+std::to_string(k + 1)+" ("+std::to_string(axes[k])+" pixels)");
                if(size[k] == 0 || origin[k] + size[k] > axes[k])
                    size[k] = axes[k] - origin[k];
            }
        }
    }

#pragma region - Description of ndarrays in their file
    namespace asdfdetail
    {
        std::shared_ptr<ArraySource> describeArray(const std::shared_ptr<FileSource>& file, asdf_value_t* value,
                                                   const std::string& path, ASDFnode desc)
        {
            auto source = std::make_shared<ArraySource>();
            source->file = file;
            source->path = path;
            source->desc = std::move(desc);

            // An integer source is a block of the file; a string one is another file
            if(const ASDFnode* block = source->desc.Find("source"); block && block->IsInteger())
                source->block = block->Get<int64_t>();

            asdf_ndarray_t* raw = nullptr;
            if(asdf_value_as_ndarray(value, &raw) == ASDF_VALUE_OK && raw)
            {
                NdarrayPtr ndarray(raw);
                source->shape.assign(raw->shape, raw->shape + raw->ndim);
                source->type = raw->datatype.type;
                const uint64_t elements = asdf_ndarray_size(raw);
                source->elementBytes = elements ? asdf_ndarray_nbytes(raw) / elements : asdf_scalar_datatype_size(raw->datatype.type);
            }
            else if(const ASDFnode* shape = source->desc.Find("shape"); shape && shape->IsSequence())
            {
                for(const ASDFnode& length : shape->Items())
                {
                    if(!length.IsInteger())
                    {
                        source->shape.clear();  // a streamed array ('*'): unknown
                        break;
                    }
                    source->shape.push_back(length.Get<uint64_t>());
                }
            }

            switch(source->type)
            {
                case ASDF_DATATYPE_STRUCTURED: source->datatype = "structured"; break;
                case ASDF_DATATYPE_ASCII:      source->datatype = "ascii";      break;
                case ASDF_DATATYPE_UCS4:       source->datatype = "ucs4";       break;
                case ASDF_DATATYPE_UNKNOWN:
                {
                    const ASDFnode* datatype = source->desc.Find("datatype");
                    source->datatype = (datatype && datatype->IsString()) ? datatype->Get<std::string>() : "unknown";
                    break;
                }
                default:
                    source->datatype = datatypeName(source->type);
                    break;
            }
            return source;
        }

        ValuePtr valueAt(asdf_file_t* file, const std::string& path)
        {
            ValuePtr current(asdf_get_value(file, ""));
            for(const std::string& part : ASDFnode::SplitPath(path))
            {
                if(!current)
                    return nullptr;

                asdf_mapping_t* mapping = nullptr;
                asdf_sequence_t* sequence = nullptr;
                if(asdf_value_as_mapping(current.get(), &mapping) == ASDF_VALUE_OK)
                    current.reset(asdf_mapping_get(mapping, part.c_str()));
                else if(asdf_value_as_sequence(current.get(), &sequence) == ASDF_VALUE_OK)
                {
                    int index = 0;
                    const auto [end, ec] = std::from_chars(part.data(), part.data() + part.size(), index);
                    if(ec != std::errc() || end != part.data() + part.size() || index < 0 || index >= asdf_sequence_size(sequence))
                        return nullptr;
                    current.reset(asdf_sequence_get(sequence, index));
                }
                else
                    return nullptr;
            }
            return current;
        }

        std::shared_ptr<ASDFarray> ArrayAccess::make(std::shared_ptr<const ArraySource> source)
        {
            std::shared_ptr<ASDFarray> array(new ASDFarray());
            array->fsource = std::move(source);
            return array;
        }

        Tile ArrayAccess::tile(const ASDFarray& array)
        {
            Tile tile;
            tile.origin  = array.forigin;
            tile.size    = array.fsize;
            tile.dropped = array.fdropped;
            return tile;
        }

        NdarrayPtr ArrayAccess::open(const ArraySource& source)
        {
            ValuePtr value = valueAt(source.file->file.get(), source.path);
            if(!value)
                throw FITSexception(READ_ERROR, kWho, "read", "no ndarray at '"+source.path+"' in "+source.file->path);

            asdf_ndarray_t* raw = nullptr;
            if(asdf_value_as_ndarray(value.get(), &raw) != ASDF_VALUE_OK || !raw)
                throw FITSexception(READ_ERROR, kWho, "read", "can't read the ndarray at '"+source.path+"' in "+source.file->path+": "+lastError(source.file->file.get()));
            return NdarrayPtr(raw);
        }
    }

#pragma endregion
#pragma region - ASDFarray
    ASDFarray::ASDFarray() = default;
    ASDFarray::~ASDFarray() = default;

#pragma region * Description
    std::vector<uint64_t> ASDFarray::Shape() const
    {
        if(fmemoryOnly || !fsource)
        {
            std::vector<uint64_t> shape;
            if(fcube)
                for(size_t i = fcube->GetDimension(); i > 0; i--)
                    shape.push_back(fcube->Size(i));
            return shape;
        }

        const std::vector<uint64_t>& full = fsize.empty() ? fsource->shape : fsize;
        if(fdropped >= full.size())
            return std::vector<uint64_t>(1, 1);
        return std::vector<uint64_t>(full.begin() + static_cast<std::ptrdiff_t>(fdropped), full.end());
    }

    std::vector<size_t> ASDFarray::Axes() const
    {
        const std::vector<uint64_t> shape = Shape();
        return std::vector<size_t>(shape.rbegin(), shape.rend());
    }

    uint64_t ASDFarray::Elements() const
    {
        const std::vector<uint64_t> shape = Shape();
        if(shape.empty())
            return 0;
        return std::accumulate(shape.begin(), shape.end(), uint64_t(1), std::multiplies<uint64_t>());
    }

    uint64_t ASDFarray::Bytes() const
    {
        if(fmemoryOnly || !fsource)
            return fcube ? Elements() * asdf_scalar_datatype_size(cubeType(*fcube)) : 0;
        return Elements() * fsource->elementBytes;
    }

    std::string ASDFarray::DataType() const
    {
        if(fmemoryOnly || !fsource)
            return fcube ? datatypeName(cubeType(*fcube)) : std::string("unknown");
        return fsource->datatype;
    }

    bool ASDFarray::IsImage() const
    {
        if(fmemoryOnly || !fsource)
            return fcube != nullptr;
        return isImageType(fsource->type);
    }

    std::string ASDFarray::Compression() const
    {
        if(!fsource || fsource->block < 0)
            return std::string();
        if(fcompressionKnown)
            return fcompression;

        std::lock_guard<std::mutex> lock(fsource->file->mtx);
        asdf_block_t* block = asdf_block_open(fsource->file->file.get(), static_cast<size_t>(fsource->block));
        if(!block)
            throw FITSexception(READ_ERROR, kWho, "Compression", "can't open block "+std::to_string(fsource->block)+" of "+fsource->file->path);
        const char* compression = asdf_block_compression(block);
        fcompression = compression ? std::string(compression) : std::string();
        asdf_block_close(block);

        // The header field is padded with zeros
        fcompression.erase(fcompression.find_last_not_of(std::string(1, '\0') + " ") + 1);
        fcompressionKnown = true;
        return fcompression;
    }

    bool ASDFarray::IsInFile() const { return fsource && !fmemoryOnly; }
    bool ASDFarray::IsLoaded() const { return fcube != nullptr; }
    bool ASDFarray::HasView()  const { return !fsize.empty(); }

#pragma endregion
#pragma region * Data
    std::shared_ptr<FITScube> ASDFarray::Load()
    {
        if(fcube)
            return fcube;
        if(!fsource)
            throw FITSexception(SHARED_NULPTR, kWho, "Load", "the handle holds no array");

        const std::vector<std::string> parts = ASDFnode::SplitPath(fsource->path);
        FITShdu header;
        if(!parts.empty())
            header.ValueForKey("EXTNAME", parts.back(), fChar);

        std::lock_guard<std::mutex> lock(fsource->file->mtx);
        NdarrayPtr ndarray = ArrayAccess::open(*fsource);
        fcube = readImage(ndarray.get(), header, ArrayAccess::tile(*this));
        return fcube;
    }

    std::shared_ptr<FITScube> ASDFarray::Read(const std::vector<size_t>& origin, const std::vector<size_t>& size) const
    {
        const std::vector<size_t> axes = Axes();
        std::vector<size_t> o(origin);
        std::vector<size_t> s(size);
        fitPart(axes, o, s, "Read");

        if(fcube)
            return cubeWindow(*fcube, o, s);
        if(!fsource)
            throw FITSexception(SHARED_NULPTR, kWho, "Read", "the handle holds no array");

        // The part, in the file and in the ASDF order: the view's kept axes are the last ones
        Tile tile = ArrayAccess::tile(*this);
        if(tile.whole())
        {
            tile.origin.assign(fsource->shape.size(), 0);
            tile.size = fsource->shape;
        }
        const size_t rank = tile.size.size();
        for(size_t k = 0; k < axes.size(); k++)
        {
            tile.origin[rank - 1 - k] += o[k];
            tile.size  [rank - 1 - k]  = s[k];
        }

        const std::vector<std::string> parts = ASDFnode::SplitPath(fsource->path);
        FITShdu header;
        if(!parts.empty())
            header.ValueForKey("EXTNAME", parts.back(), fChar);

        std::lock_guard<std::mutex> lock(fsource->file->mtx);
        NdarrayPtr ndarray = ArrayAccess::open(*fsource);
        return readImage(ndarray.get(), header, tile);
    }

    double ASDFarray::ValueAt(const std::vector<size_t>& pixel) const
    {
        const std::vector<size_t> axes = Axes();
        if(pixel.size() != axes.size())
            throw FITSexception(BAD_DIMEN, kWho, "ValueAt", "the array has "+std::to_string(axes.size())+" axes, the position "+std::to_string(pixel.size()));
        for(size_t k = 0; k < axes.size(); k++)
            if(pixel[k] >= axes[k])
                throw FITSexception(BAD_DIMEN, kWho, "ValueAt", "position "+std::to_string(pixel[k])+" is past the end of axis "+std::to_string(k + 1));

        if(fcube)
        {
            size_t index = 0;
            size_t stride = 1;
            for(size_t k = 0; k < axes.size(); k++)
            {
                index += pixel[k] * stride;
                stride *= axes[k];
            }
            return fcube->at(index);
        }
        return Read(pixel, std::vector<size_t>(pixel.size(), 1))->at(0);
    }

    bool ASDFarray::Release()
    {
        if(fmemoryOnly)
            return false;
        fcube.reset();
        return true;
    }

    void ASDFarray::Set(std::shared_ptr<FITScube> cube)
    {
        if(!cube)
            throw FITSexception(SHARED_NULPTR, kWho, "Set", "no image given");
        fcube = std::move(cube);
        fmemoryOnly = true;
        forigin.clear();
        fsize.clear();
        fdropped = 0;
    }

#pragma endregion
#pragma region * Views and copies
    std::shared_ptr<ASDFarray> ASDFarray::Window(const std::vector<size_t>& origin, const std::vector<size_t>& size) const
    {
        const std::vector<size_t> axes = Axes();
        std::vector<size_t> o(origin);
        std::vector<size_t> s(size);
        fitPart(axes, o, s, "Window");

        std::shared_ptr<ASDFarray> out(new ASDFarray());
        out->fsource = fsource;
        out->fcompression = fcompression;
        out->fcompressionKnown = fcompressionKnown;
        out->fmemoryOnly = fmemoryOnly;
        if(fcube)
            out->fcube = cubeWindow(*fcube, o, s);
        if(fmemoryOnly || !fsource)
            return out;

        Tile tile = ArrayAccess::tile(*this);
        if(tile.whole())
        {
            tile.origin.assign(fsource->shape.size(), 0);
            tile.size = fsource->shape;
        }
        const size_t rank = tile.size.size();
        for(size_t k = 0; k < axes.size(); k++)
        {
            tile.origin[rank - 1 - k] += o[k];
            tile.size  [rank - 1 - k]  = s[k];
        }
        out->forigin  = std::move(tile.origin);
        out->fsize    = std::move(tile.size);
        out->fdropped = tile.dropped;
        return out;
    }

    std::shared_ptr<ASDFarray> ASDFarray::Layer(size_t index) const
    {
        const std::vector<size_t> axes = Axes();
        if(axes.size() < 2)
            throw FITSexception(BAD_DIMEN, kWho, "Layer", "the array has "+std::to_string(axes.size())+" axis: no layers");
        if(index >= axes.back())
            throw FITSexception(BAD_DIMEN, kWho, "Layer", "the array has "+std::to_string(axes.back())+" layers");

        std::vector<size_t> origin(axes.size(), 0);
        std::vector<size_t> size(axes);
        origin.back() = index;
        size.back() = 1;

        std::shared_ptr<ASDFarray> out(new ASDFarray());
        out->fsource = fsource;
        out->fcompression = fcompression;
        out->fcompressionKnown = fcompressionKnown;
        out->fmemoryOnly = fmemoryOnly;
        if(fcube)
            out->fcube = cubeWindow(*fcube, origin, size, 1);
        if(fmemoryOnly || !fsource)
            return out;

        Tile tile = ArrayAccess::tile(*this);
        if(tile.whole())
        {
            tile.origin.assign(fsource->shape.size(), 0);
            tile.size = fsource->shape;
        }
        // The last FITS axis is the first one the view keeps
        tile.origin[tile.dropped] += index;
        tile.size  [tile.dropped]  = 1;
        out->forigin  = std::move(tile.origin);
        out->fsize    = std::move(tile.size);
        out->fdropped = tile.dropped + 1;
        return out;
    }

    std::shared_ptr<ASDFarray> ASDFarray::Clone() const
    {
        std::shared_ptr<ASDFarray> out(new ASDFarray());
        out->fsource = fsource;
        out->forigin = forigin;
        out->fsize = fsize;
        out->fdropped = fdropped;
        out->fmemoryOnly = fmemoryOnly;
        out->fcompression = fcompression;
        out->fcompressionKnown = fcompressionKnown;
        if(fcube)
            out->fcube = copyCube(*fcube);
        return out;
    }

    std::shared_ptr<ASDFarray> ASDFarray::FromCube(std::shared_ptr<FITScube> cube)
    {
        if(!cube)
            throw FITSexception(SHARED_NULPTR, kWho, "FromCube", "no image given");
        std::shared_ptr<ASDFarray> out(new ASDFarray());
        out->fcube = std::move(cube);
        out->fmemoryOnly = true;
        return out;
    }
#pragma endregion
#pragma endregion
}
