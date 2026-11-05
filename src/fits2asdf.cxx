//
//  fits2asdf.cxx
//  asdf2fits
//
//  Written by Claude (Anthropic's AI model) for William Gillard.
//  Copyright (c) 2026 William Gillard.
//
//  SPDX-License-Identifier: BSD-3-Clause
//  Licence text: see the LICENSE file at the root of this project.

/**
 * @file fits2asdf.cxx
 * @brief Implementation of fits2asdf(): an extension of libasdf serializing a list of HDUs.
 */
#include "fits2asdf.h"

#include <exception>
#include <filesystem>
#include <mutex>
#include <vector>

#include "ASDFheader.h"
#include "ASDFhandle.h"
#include "ASDFimage.h"
#include "ASDFtable.h"

namespace DSL
{
    using namespace asdfdetail;

    namespace
    {
        const char* const kWho = "fits2asdf";   //!< Name given in the exceptions thrown here

        /** @brief One HDU to write. With neither image nor table, only its header is written. */
        struct HDUOut
        {
            const FITScube*  cube  = nullptr;   //!< Image to write, not owned; may be null
            const FITStable* table = nullptr;   //!< Table to write, not owned; may be null
            FITShdu header;                     //!< Copy of the header, adjusted while writing
        };

        /**
         * @brief What the extension serializes.
         * @details Exceptions can't cross the C frames of libasdf: they are kept here and thrown
         * again once back in C++.
         */
        struct HDUList
        {
            std::vector<HDUOut> hdus;       //!< HDUs in file order; the first is the primary HDU
            std::exception_ptr failure;     //!< Exception caught while serializing, if any
        };

        /** @brief Structural keywords, written first in this order as in a FITS header. */
        const std::vector<std::string> kFirstCards = {"SIMPLE", "XTENSION", "BITPIX", "NAXIS", "NAXIS1", "NAXIS2", "NAXIS3", "NAXIS4",
                                                      "EXTEND", "PCOUNT", "GCOUNT", "TFIELDS"};

        /**
         * @brief Remove the keywords that are set again on writing or no longer hold.
         * @param header Header to update: SIMPLE, XTENSION, EXTEND, PCOUNT, GCOUNT, CHECKSUM and DATASUM go.
         */
        void dropStaleKeywords(FITShdu& header)
        {
            for(const char* key : {"SIMPLE", "XTENSION", "EXTEND", "PCOUNT", "GCOUNT", "CHECKSUM", "DATASUM"})
                if(header.Exists(key))
                    header.DeleteKey(key);
        }

        /**
         * @brief Turn one HDU into a {header, data} mapping.
         * @param file    File the mapping is created in.
         * @param hdu     HDU to write; its header is adjusted to match what is written.
         * @param primary True for the primary HDU.
         * @return The new mapping; it has no @c data entry if the HDU holds no data.
         * @throws FITSexception As writeImage() or writeTable().
         */
        asdf_value_t* serializeHDU(asdf_file_t* file, HDUOut& hdu, bool primary)
        {
            asdf_mapping_t* node = asdf_mapping_create(file);

            asdf_value_t* data = nullptr;
            dropStaleKeywords(hdu.header);
            if(hdu.cube)
                data = writeImage(file, *hdu.cube, hdu.header);
            else if(hdu.table)
                data = writeTable(file, *hdu.table, hdu.header);
            else
            {
                clearImageStructure(hdu.header);
                hdu.header.ValueForKey("BITPIX", int32_t(8));
                hdu.header.ValueForKey("NAXIS", int32_t(0));
            }

            if(primary)
            {
                hdu.header.ValueForKey("SIMPLE", true, "conforms to FITS standard");
                hdu.header.ValueForKey("EXTEND", true);
            }
            else if(!hdu.table)
            {
                hdu.header.ValueForKey("XTENSION", std::string("IMAGE"), fChar, "image extension");
                hdu.header.ValueForKey("PCOUNT", uint64_t(0));
                hdu.header.ValueForKey("GCOUNT", uint64_t(1));
            }

            if(data)
                asdf_mapping_set(node, "data", data);

            asdf_mapping_set_sequence(node, "header", writeCards(file, hdu.header, kFirstCards));
            return asdf_value_of_mapping(node);
        }

        /**
         * @brief Serialize callback of the extension: turn the HDU list into a sequence.
         * @param file   File the sequence is created in.
         * @param object The HDUList to write. On an exception, its @c failure member keeps it.
         * @return The new sequence, or null if an exception was caught.
         */
        asdf_value_t* serializeList(asdf_file_t* file, const void* object, const void*)
        {
            auto* list = const_cast<HDUList*>(static_cast<const HDUList*>(object));
            try
            {
                asdf_sequence_t* sequence = asdf_sequence_create(file);
                for(size_t i = 0; i < list->hdus.size(); i++)
                    asdf_sequence_append(sequence, serializeHDU(file, list->hdus[i], i == 0));
                return asdf_value_of_sequence(sequence);
            }
            catch(...)
            {
                list->failure = std::current_exception();
                return nullptr;
            }
        }

        /**
         * @brief Get the libasdf extension for the astropy FITS tag, registering it on first use.
         * @return The extension, valid until the end of the program.
         */
        asdf_extension_t* fitsExtension()
        {
            static const char* const tags[] = {"tag:astropy.org:astropy/fits/fits-1.2.0", nullptr};
            static asdf_extension_vtab_t vtab{};
            static asdf_extension_t extension{};
            static std::once_flag once;

            std::call_once(once, []{
                vtab.serialize = serializeList;
                extension.tags = tags;
                extension.vtab = &vtab;
                extension.size = sizeof(HDUList);
                asdf_extension_register(&extension);
            });
            return &extension;
        }

        /**
         * @brief Serialize an HDU list and write the file.
         * @param list    HDUs to write.
         * @param path    File to create.
         * @param replace Overwrite an existing file.
         * @throws FITSexception FILE_NOT_CREATED if the file exists and @p replace is false, or if it can't
         * be created or written. WRITE_ERROR if serializing fails, or the exception caught while serializing.
         */
        void writeFile(HDUList& list, const std::string& path, bool replace)
        {
            if(!replace && std::filesystem::exists(path))
                throw FITSexception(FILE_NOT_CREATED, kWho, "write", path+" already exists");

            FilePtr file(asdf_open_mem(nullptr, 0));
            if(!file)
                throw FITSexception(FILE_NOT_CREATED, kWho, "write", "can't create an ASDF file: "+lastError(nullptr));

            const asdf_value_err_t err = asdf_set_extension_type(file.get(), "fits", &list, fitsExtension());
            if(list.failure)
                std::rethrow_exception(list.failure);
            if(err != ASDF_VALUE_OK)
                throw FITSexception(WRITE_ERROR, kWho, "write", lastError(file.get()));

            if(asdf_write_to_file(file.get(), path.c_str()) != 0)
                throw FITSexception(FILE_NOT_CREATED, kWho, "write", "can't write "+path+": "+lastError(file.get()));
        }

        /**
         * @brief Describe an image to write.
         * @param cube Image, which must outlive the result.
         * @return The HDU, with a copy of the image's header.
         */
        HDUOut cubeHDU(const FITScube& cube)
        {
            HDUOut out;
            out.cube = &cube;
            out.header = cube.HDU();
            return out;
        }

        /**
         * @brief Describe a table to write.
         * @param table Table, which must outlive the result.
         * @return The HDU, with a copy of the table's header.
         */
        HDUOut tableHDU(const FITStable& table)
        {
            HDUOut out;
            out.table = &table;
            out.header = table.HDU();
            return out;
        }
    }

    void fits2asdf(const FITScube& cube, const std::string& path, bool replace)
    {
        HDUList list;
        list.hdus.push_back(cubeHDU(cube));
        writeFile(list, path, replace);
    }

    void fits2asdf(const FITStable& table, const std::string& path, bool replace)
    {
        HDUList list;
        list.hdus.emplace_back();
        list.hdus.push_back(tableHDU(table));
        writeFile(list, path, replace);
    }

    void fits2asdf(FITSmanager& fits, const std::string& path, bool replace)
    {
        // The objects have to outlive the list that refers to them
        std::vector<std::shared_ptr<FITScube>>  cubes;
        std::vector<std::shared_ptr<FITStable>> tables;

        HDUList list;
        for(int i = 1; i <= fits.NumberOfHeader(); i++)
        {
            const int type = fits.MoveToHDU(i);
            if(type == IMAGE_HDU)
            {
                std::shared_ptr<FITShdu> header = fits.GetHeaderAtIndex(i);
                if(header->GetDimension() > 0)
                {
                    cubes.push_back(fits.GetImageAtIndex(i));
                    list.hdus.push_back(cubeHDU(*cubes.back()));
                }
                else
                {
                    HDUOut empty;
                    empty.header = *header;
                    list.hdus.push_back(empty);
                }
            }
            else
            {
                tables.push_back(fits.GetTableAtIndex(i));
                list.hdus.push_back(tableHDU(*tables.back()));
            }
        }

        writeFile(list, path, replace);
    }
}
