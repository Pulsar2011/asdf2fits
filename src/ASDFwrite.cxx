//
//  ASDFwrite.cxx
//  asdf2fits
//
//  Written by Claude (Anthropic's AI model) for William Gillard.
//  Copyright (c) 2026 William Gillard.
//
//  SPDX-License-Identifier: BSD-3-Clause
//  Licence text: see the LICENSE file at the root of this project.

/**
 * @file ASDFwrite.cxx
 * @brief Implementation of the writing of an ASDFnode tree, with its arrays, as an ASDF file.
 * @details The tree is written as YAML text, which libasdf parses into a new file; the blocks are then
 * appended in the order the YAML refers to them, and libasdf writes the file. Blocks copied from another file
 * borrow its memory-mapped bytes, so they are neither decompressed nor copied in memory.
 */
#include "ASDFwrite.h"

#include <algorithm>
#include <bit>
#include <chrono>
#include <ctime>
#include <filesystem>
#include <map>
#include <mutex>
#include <random>
#include <set>

#include "ASDFsource.h"
#include "ASDFtree.h"

#ifndef ASDF2FITS_VERSION
#define ASDF2FITS_VERSION "unknown"     //!< Version written in history entries; set by the build
#endif

namespace DSL::asdfdetail
{
    namespace
    {
        const char* const kWho = "WriteASDF";   //!< Name given in the exceptions thrown here

        /** @brief One block of the file to write. */
        struct BlockOut
        {
            PixelBlock pixels;                  //!< Pixels written from memory; unused for a copy
            std::shared_ptr<FITScube> keep;     //!< Image owning the pixels, kept alive until the file is written
            std::string compression;            //!< Compression of pixels written from memory
            std::shared_ptr<FileSource> file;   //!< File of a block copied as it is; null for pixels from memory
            int64_t index = -1;                 //!< Index of the copied block in @c file
        };

        /** @brief Closes the blocks opened in other files, once the file is written. */
        struct OpenedBlocks
        {
            std::vector<asdf_block_t*> blocks;  //!< Blocks to close
            ~OpenedBlocks() { for(asdf_block_t* block : blocks) asdf_block_close(block); }
        };

        /** @brief Removes the temporary file unless it was renamed. */
        struct TemporaryFile
        {
            std::filesystem::path path;         //!< The temporary file
            bool kept = false;                  //!< Renamed to the destination
            ~TemporaryFile()
            {
                std::error_code ec;
                if(!kept)
                    std::filesystem::remove(path, ec);
            }
        };

        /**
         * @brief Trim a compression name read from a block header.
         * @param name The name, possibly padded.
         * @return The name, or empty for none.
         */
        std::string compressionName(const char* name)
        {
            std::string text = name ? std::string(name) : std::string();
            text.erase(text.find_last_not_of(std::string(1, '\0') + " ") + 1);
            return text;
        }

        /**
         * @brief Make the description of an ndarray written from memory.
         * @param pixels Its pixels.
         * @param index  Index of its block.
         * @return The tagged mapping; inline and empty if the array has no element.
         */
        ASDFnode ndarrayOf(const PixelBlock& pixels, int64_t index)
        {
            ASDFnode desc = ASDFnode::Mapping(kNdarrayTag);
            if(pixels.bytes == 0)
                desc["data"] = ASDFnode::Sequence();
            else
                desc["source"] = index;
            desc["datatype"] = datatypeName(pixels.type);
            desc["byteorder"] = std::string(std::endian::native == std::endian::little ? "little" : "big");
            ASDFnode& shape = desc["shape"];
            shape = ASDFnode::Sequence();
            for(uint64_t length : pixels.shape)
                shape.Append(length);
            return desc;
        }

        /**
         * @brief Get the time now, as written in history entries.
         * @return The UTC time, as "YYYY-MM-DDThh:mm:ss".
         */
        std::string utcNow()
        {
            const std::time_t now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
            std::tm utc{};
            gmtime_r(&now, &utc);
            char text[32];
            std::strftime(text, sizeof(text), "%Y-%m-%dT%H:%M:%S", &utc);
            return text;
        }

        /**
         * @brief Make a name for the temporary file next to the destination.
         * @param path The destination.
         * @return A hidden name in the same directory, unlikely to exist.
         */
        std::filesystem::path temporaryName(const std::filesystem::path& path)
        {
            std::random_device device;
            const std::string suffix = std::to_string(device()) + std::to_string(device());
            return path.parent_path() / ("." + path.filename().string() + ".tmp" + suffix);
        }
    }

    void checkCompression(const std::string& compression)
    {
        static const std::set<std::string> known = {"keep", "", "zlib", "bzp2", "lz4"};
        if(!known.count(compression))
            throw FITSexception(BAD_OPTION, kWho, "options", "unknown compression '"+compression+"' (keep, zlib, bzp2, lz4 or none)");
    }

    ASDFnode historyEntry(const std::string& description)
    {
        ASDFnode software = ASDFnode::Mapping("tag:stsci.edu:asdf/core/software-1.0.0");
        software["name"] = "asdf2fits";
        software["version"] = ASDF2FITS_VERSION;
        software["homepage"] = "https://gitlab.in2p3.fr/deepskyproject/asdf2fits";

        ASDFnode entry = ASDFnode::Mapping("tag:stsci.edu:asdf/core/history_entry-1.0.0");
        entry["description"] = description;
        entry["time"] = utcNow();
        entry["software"] = std::move(software);
        return entry;
    }

    void appendHistory(ASDFnode& tree, ASDFnode entry)
    {
        ASDFnode& history = tree["history"];
        if(history.IsSequence())
        {
            history.Append(std::move(entry));   // the old style: a list of entries
            return;
        }
        if(!history.IsMapping())
            history = ASDFnode::Mapping();
        ASDFnode& entries = history["entries"];
        if(!entries.IsSequence())
            entries = ASDFnode::Sequence();
        entries.Append(std::move(entry));
    }

    void writeTreeFile(const ASDFnode& tree, const std::string& imagePath, PixelBlock image,
                       const std::string& compression, const std::string& path, bool replace)
    {
        checkCompression(compression);
        if(!replace && std::filesystem::exists(path))
            throw FITSexception(FILE_NOT_CREATED, kWho, "write", path+" already exists");

        // The YAML, and what each block it refers to holds
        std::vector<BlockOut> blocks;
        std::map<std::pair<const FileSource*, int64_t>, int64_t> copies;

        auto fromMemory = [&](PixelBlock pixels, std::shared_ptr<FITScube> keep, const std::string& original)
        {
            BlockOut block;
            block.pixels = std::move(pixels);
            block.keep = std::move(keep);
            block.compression = (compression == "keep") ? original : compression;
            const int64_t index = static_cast<int64_t>(blocks.size());
            ASDFnode desc = ndarrayOf(block.pixels, index);
            if(block.pixels.bytes > 0)
                blocks.push_back(std::move(block));
            return desc;
        };

        const ArrayWriter arrays = [&](const std::string& nodePath, const ASDFnode& node) -> ASDFnode
        {
            const ASDFarray& array = *node.Array();
            const ArraySource* source = ArrayAccess::source(array);
            const std::string original = source ? array.Compression() : std::string();

            if(!imagePath.empty() && nodePath == imagePath)
                return fromMemory(std::move(image), nullptr, original);

            if(const std::shared_ptr<FITScube>& cube = ArrayAccess::cube(array))
                return fromMemory(cubePixels(*cube), cube, original);

            if(!source)
                return ASDFnode();

            if(array.HasView())
            {
                std::shared_ptr<FITScube> tile = array.Read({}, {});
                return fromMemory(cubePixels(*tile), tile, original);
            }

            ASDFnode desc = source->desc;
            if(source->block >= 0)
            {
                // A block copied as it is; arrays sharing a block keep sharing it
                const auto key = std::make_pair(source->file.get(), source->block);
                auto found = copies.find(key);
                if(found == copies.end())
                {
                    BlockOut block;
                    block.file = source->file;
                    block.index = source->block;
                    found = copies.emplace(key, static_cast<int64_t>(blocks.size())).first;
                    blocks.push_back(std::move(block));
                }
                desc["source"] = found->second;
            }
            return desc;    // inline and external arrays stay as they were
        };

        const std::string yaml = writeDocument(tree, arrays);

        FilePtr file(asdf_open_mem(yaml.data(), yaml.size()));
        if(!file || asdf_error(file.get()))
            throw FITSexception(WRITE_ERROR, kWho, "write", "libasdf refuses the tree: "+lastError(file.get()));

        // The files blocks are copied from stay locked until the file is written: lock them in address order
        std::set<FileSource*> sources;
        for(const BlockOut& block : blocks)
            if(block.file)
                sources.insert(block.file.get());
        std::vector<std::unique_lock<std::mutex>> locks;
        for(FileSource* source : sources)
            locks.emplace_back(source->mtx);

        OpenedBlocks opened;
        for(const BlockOut& block : blocks)
        {
            asdf_block_t* out = nullptr;
            if(block.file)
            {
                asdf_block_t* in = asdf_block_open(block.file->file.get(), static_cast<size_t>(block.index));
                if(!in)
                    throw FITSexception(READ_ERROR, kWho, "write", "can't open block "+std::to_string(block.index)+" of "+block.file->path);
                opened.blocks.push_back(in);

                size_t size = 0;
                const void* raw = asdf_block_data_raw(in, &size);
                if(!raw && size > 0)
                    throw FITSexception(READ_ERROR, kWho, "write", "can't read block "+std::to_string(block.index)+" of "+block.file->path);

                const std::string inCompression = compressionName(asdf_block_compression(in));
                if(inCompression.empty())
                    out = asdf_block_create(raw, size);
                else
                {
                    out = asdf_block_create(nullptr, 0);
                    if(out && asdf_block_data_set_compressed(out, raw, size, asdf_block_data_size(in), inCompression.c_str()) != 0)
                    {
                        asdf_block_destroy(out);
                        throw FITSexception(WRITE_ERROR, kWho, "write", "can't copy block "+std::to_string(block.index)+" of "+block.file->path);
                    }
                }
            }
            else
            {
                out = asdf_block_create(block.pixels.data, block.pixels.bytes);
                if(out && !block.compression.empty() && asdf_block_compression_set(out, block.compression.c_str()) != 0)
                {
                    asdf_block_destroy(out);
                    throw FITSexception(WRITE_ERROR, kWho, "write", "libasdf can't compress with "+block.compression);
                }
            }

            if(!out)
                throw FITSexception(MEMORY_ALLOCATION, kWho, "write", "can't create a block: "+lastError(file.get()));
            if(!asdf_block_append(file.get(), out))
            {
                asdf_block_destroy(out);
                throw FITSexception(WRITE_ERROR, kWho, "write", "can't add a block: "+lastError(file.get()));
            }
        }

        writeReplacing(file.get(), path);
    }

    void writeReplacing(asdf_file_t* file, const std::string& path)
    {
        const std::filesystem::path destination(path);
        TemporaryFile temporary{temporaryName(destination)};
        if(asdf_write_to_file(file, temporary.path.c_str()) != 0)
            throw FITSexception(FILE_NOT_CREATED, kWho, "write", "can't write "+path+": "+lastError(file));

        std::error_code ec;
        std::filesystem::rename(temporary.path, destination, ec);
        if(ec)
            throw FITSexception(FILE_NOT_CREATED, kWho, "write", "can't rename the file written to "+path+": "+ec.message());
        temporary.kept = true;
    }
}
