//
//  ASDFmanager.cxx
//  asdf2fits
//
//  Written by Claude (Anthropic's AI model) for William Gillard.
//  Copyright (c) 2026 William Gillard.
//
//  SPDX-License-Identifier: BSD-3-Clause
//  Licence text: see the LICENSE file at the root of this project.

#include "ASDFmanager.h"

#include <algorithm>
#include <cctype>
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
        const char* const kWho = "ASDFmanager";

        std::string lower(std::string text)
        {
            std::transform(text.begin(), text.end(), text.begin(), [](unsigned char c){ return std::tolower(c); });
            return text;
        }

        struct HDUEntry
        {
            enum Kind { empty, image, table } kind = empty;
            std::string name;
            FITShdu header;
            NdarrayPtr ndarray;                 //!< Image pixels, or the rows of a structured table
            std::vector<NdarrayPtr> arrays;     //!< One per column of a column-based table
            std::vector<ColumnSpec> columns;
            uint64_t rows = 0;
            bool structured = false;
        };

        bool isTableValue(asdf_value_t* value)
        {
            return hasTag(value, "table/table-") || hasTag(value, "core/table-");
        }

        // True when the value holds, at any depth, something that becomes an HDU
        bool holdsData(asdf_value_t* value)
        {
            if(isOpaqueValue(value))
                return false;
            if(isDataValue(value))
                return true;

            asdf_container_iter_t* it = asdf_value_is_container(value) ? asdf_container_iter_init(value) : nullptr;
            while(it && asdf_container_iter_next(&it))
            {
                if(holdsData(it->value))
                {
                    asdf_container_iter_destroy(it);
                    return true;
                }
            }
            return false;
        }

        void describeData(asdf_value_t* value, HDUEntry& entry)
        {
            if(isTableValue(value))
            {
                entry.kind = HDUEntry::table;
                entry.columns = columnTable(value, entry.arrays);
                entry.rows = entry.arrays.empty() ? 0 : entry.arrays.front()->shape[0];
                for(const NdarrayPtr& array : entry.arrays)
                    if(array->shape[0] != entry.rows)
                        throw FITSexception(BAD_DIMEN, kWho, "open", "the columns of table '"+entry.name+"' have different lengths");

                setTableStructure(entry.header, entry.columns, entry.rows);
                return;
            }

            asdf_ndarray_t* raw = nullptr;
            if(asdf_value_as_ndarray(value, &raw) != ASDF_VALUE_OK)
                throw FITSexception(BAD_DATATYPE, kWho, "open", "can't read the ndarray of HDU '"+entry.name+"'");
            entry.ndarray.reset(raw);

            if(raw->datatype.type == ASDF_DATATYPE_STRUCTURED)
            {
                entry.kind = HDUEntry::table;
                entry.structured = true;
                entry.columns = structuredColumns(raw);
                entry.rows = raw->ndim > 0 ? raw->shape[0] : 0;
                setTableStructure(entry.header, entry.columns, entry.rows);
                return;
            }

            if(asdf_ndarray_size(raw) > 0 && raw->ndim > 0)
            {
                entry.kind = HDUEntry::image;
                setImageStructure(entry.header, raw->datatype.type, fitsAxes(raw));
            }
            else
            {
                clearImageStructure(entry.header);
                entry.header.ValueForKey("BITPIX", int32_t(8));
                entry.header.ValueForKey("NAXIS", int32_t(0));
            }
        }
    }

    struct ASDFmanager::Impl
    {
        mutable std::mutex mtx;
        FilePtr file;
        std::string path;
        std::vector<HDUEntry> hdus;
        bool fitsLayout = false;
        FITShdu rootKeywords;

        void addFitsLayout(asdf_value_t* value)
        {
            fitsLayout = true;
            asdf_sequence_t* list = nullptr;
            if(asdf_value_as_sequence(value, &list) != ASDF_VALUE_OK)
                return;

            asdf_sequence_iter_t* it = asdf_sequence_iter_init(list);
            while(asdf_sequence_iter_next(&it))
            {
                asdf_mapping_t* node = nullptr;
                if(asdf_value_as_mapping(it->value, &node) != ASDF_VALUE_OK)
                    continue;

                HDUEntry entry;
                ValuePtr cards(asdf_mapping_get(node, "header"));
                if(cards)
                    readCards(cards.get(), entry.header);

                if(entry.header.Exists("EXTNAME"))
                    entry.name = entry.header.GetValueForKey("EXTNAME");
                else
                    entry.name = hdus.empty() ? "PRIMARY" : "HDU" + std::to_string(hdus.size() + 1);

                ValuePtr data(asdf_mapping_get(node, "data"));
                if(data && !asdf_value_is_null(data.get()))
                    describeData(data.get(), entry);

                hdus.push_back(std::move(entry));
            }
        }

        void addData(asdf_value_t* value, const std::string& path, asdf_value_t* siblings)
        {
            HDUEntry entry;
            entry.name = path;
            if(siblings)
                readTreeKeywords(siblings, std::string(), entry.header);
            if(lower(path) != "primary")
                entry.header.ValueForKey("EXTNAME", path, fChar);
            describeData(value, entry);
            hdus.push_back(std::move(entry));
        }

        void walk(asdf_value_t* value, const std::string& path, asdf_value_t* siblings)
        {
            if(isOpaqueValue(value))
                return;
            if(hasTag(value, "/fits/fits-"))
            {
                addFitsLayout(value);
                return;
            }
            if(isDataValue(value))
            {
                addData(value, path, siblings);
                return;
            }

            const std::string prefix = path.empty() ? std::string() : path + ".";
            asdf_mapping_t* mapping = nullptr;
            asdf_sequence_t* sequence = nullptr;
            if(asdf_value_as_mapping(value, &mapping) == ASDF_VALUE_OK)
            {
                asdf_mapping_iter_t* it = asdf_mapping_iter_init(mapping);
                while(asdf_mapping_iter_next(&it))
                {
                    if(path.empty() && (std::string(it->key) == "asdf_library" || std::string(it->key) == "history"))
                        continue;
                    if(holdsData(it->value))
                        walk(it->value, prefix + it->key, path.empty() ? nullptr : value);
                }
            }
            else if(asdf_value_as_sequence(value, &sequence) == ASDF_VALUE_OK)
            {
                asdf_sequence_iter_t* it = asdf_sequence_iter_init(sequence);
                while(asdf_sequence_iter_next(&it))
                    if(holdsData(it->value))
                        walk(it->value, prefix + std::to_string(it->index), siblings);
            }
        }

        void collectRootKeywords(asdf_value_t* root)
        {
            asdf_mapping_t* mapping = nullptr;
            if(asdf_value_as_mapping(root, &mapping) != ASDF_VALUE_OK)
                return;

            asdf_mapping_iter_t* it = asdf_mapping_iter_init(mapping);
            while(asdf_mapping_iter_next(&it))
            {
                const std::string key(it->key);
                if(key == "asdf_library" || key == "history" || holdsData(it->value) || isOpaqueValue(it->value))
                    continue;

                std::string upper(key);
                std::transform(upper.begin(), upper.end(), upper.begin(), [](unsigned char c){ return std::toupper(c); });
                if(asdf_value_is_mapping(it->value))
                    readTreeKeywords(it->value, upper, rootKeywords);
                else
                    readScalar(it->value, upper, rootKeywords);
            }
        }

        void discover()
        {
            ValuePtr root(asdf_get_value(file.get(), ""));
            if(!root)
                throw FITSexception(FILE_NOT_OPENED, kWho, "open", "no ASDF tree in "+path+": "+lastError(file.get()));

            walk(root.get(), std::string(), nullptr);

            // The primary HDU comes first, as in a FITS file
            auto primary = std::find_if(hdus.begin(), hdus.end(), [](const HDUEntry& e){ return lower(e.name) == "primary"; });
            if(primary != hdus.end())
                std::rotate(hdus.begin(), primary, primary + 1);

            if(fitsLayout)
                return;

            if(hdus.empty())
            {
                HDUEntry entry;
                entry.name = "PRIMARY";
                entry.header.ValueForKey("BITPIX", int32_t(8));
                entry.header.ValueForKey("NAXIS", int32_t(0));
                hdus.push_back(std::move(entry));
            }

            collectRootKeywords(root.get());
            hdus.front().header.Import(rootKeywords);
        }

        void requireOpen(const char* what) const
        {
            if(!file)
                throw FITSexception(FILE_NOT_OPENED, kWho, what, "the file is closed");
        }

        HDUEntry& at(int index, const char* what)
        {
            requireOpen(what);
            if(index < 1 || index > static_cast<int>(hdus.size()))
                throw FITSexception(BAD_HDU_NUM, kWho, what, "HDU "+std::to_string(index)+" doesn't exist in "+path+" ("+std::to_string(hdus.size())+" HDU)");
            return hdus[static_cast<size_t>(index - 1)];
        }

        HDUEntry& named(const std::string& name, const char* what)
        {
            requireOpen(what);
            auto it = std::find_if(hdus.begin(), hdus.end(), [&](const HDUEntry& e){ return lower(e.name) == lower(name); });
            if(it == hdus.end())
                throw FITSexception(BAD_HDU_NUM, kWho, what, "no HDU named '"+name+"' in "+path);
            return *it;
        }

        std::shared_ptr<FITScube> loadImage(HDUEntry& entry)
        {
            if(entry.kind == HDUEntry::table)
                throw FITSexception(NOT_IMAGE, kWho, "GetImage", "HDU '"+entry.name+"' of "+path+" is a table");
            if(entry.kind == HDUEntry::empty)
                throw FITSexception(NOT_IMAGE, kWho, "GetImage", "HDU '"+entry.name+"' of "+path+" holds no data");
            return readImage(entry.ndarray.get(), entry.header);
        }

        std::shared_ptr<FITStable> loadTable(HDUEntry& entry)
        {
            if(entry.kind != HDUEntry::table)
                throw FITSexception(NOT_TABLE, kWho, "GetTable", "HDU '"+entry.name+"' of "+path+" isn't a table");
            return readTable(entry.columns, entry.rows, entry.structured, entry.header);
        }
    };

    ASDFmanager::ASDFmanager(const std::string& path): impl(std::make_unique<Impl>())
    {
        impl->path = path;
        impl->file.reset(asdf_open_file(path.c_str(), "r"));
        if(!impl->file)
            throw FITSexception(FILE_NOT_OPENED, kWho, "ctor", "can't open "+path+": "+lastError(nullptr));
        if(asdf_error(impl->file.get()))
            throw FITSexception(FILE_NOT_OPENED, kWho, "ctor", "can't open "+path+": "+lastError(impl->file.get()));

        impl->discover();
    }

    ASDFmanager::~ASDFmanager() = default;

    int ASDFmanager::NumberOfHeader() const
    {
        std::lock_guard<std::mutex> lk(impl->mtx);
        return static_cast<int>(impl->hdus.size());
    }

    bool ASDFmanager::isOpen() const
    {
        std::lock_guard<std::mutex> lk(impl->mtx);
        return static_cast<bool>(impl->file);
    }

    void ASDFmanager::Close()
    {
        std::lock_guard<std::mutex> lk(impl->mtx);
        // The ndarrays depend on the file: they go first
        impl->hdus.clear();
        impl->file.reset();
    }

    std::shared_ptr<FITShdu> ASDFmanager::GetPrimaryHeader()
    {
        return GetHeaderAtIndex(1);
    }

    std::shared_ptr<FITShdu> ASDFmanager::GetHeaderAtIndex(const int& index)
    {
        std::lock_guard<std::mutex> lk(impl->mtx);
        return std::make_shared<FITShdu>(impl->at(index, "GetHeaderAtIndex").header);
    }

    std::shared_ptr<FITShdu> ASDFmanager::GetHeader(const std::string& name)
    {
        std::lock_guard<std::mutex> lk(impl->mtx);
        return std::make_shared<FITShdu>(impl->named(name, "GetHeader").header);
    }

    std::shared_ptr<FITScube> ASDFmanager::GetPrimary()
    {
        return GetImageAtIndex(1);
    }

    std::shared_ptr<FITScube> ASDFmanager::GetImageAtIndex(const int& index)
    {
        std::lock_guard<std::mutex> lk(impl->mtx);
        return impl->loadImage(impl->at(index, "GetImageAtIndex"));
    }

    std::shared_ptr<FITScube> ASDFmanager::GetImage(const std::string& name)
    {
        std::lock_guard<std::mutex> lk(impl->mtx);
        return impl->loadImage(impl->named(name, "GetImage"));
    }

    std::shared_ptr<FITStable> ASDFmanager::GetTableAtIndex(const int& index)
    {
        std::lock_guard<std::mutex> lk(impl->mtx);
        return impl->loadTable(impl->at(index, "GetTableAtIndex"));
    }

    std::shared_ptr<FITStable> ASDFmanager::GetTable(const std::string& name)
    {
        std::lock_guard<std::mutex> lk(impl->mtx);
        return impl->loadTable(impl->named(name, "GetTable"));
    }
}
