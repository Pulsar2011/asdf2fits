//
//  ASDFmanager.cxx
//  asdf2fits
//
//  Written by Claude (Anthropic's AI model) for William Gillard.
//  Copyright (c) 2026 William Gillard.
//
//  SPDX-License-Identifier: BSD-3-Clause
//  Licence text: see the LICENSE file at the root of this project.

/**
 * @file ASDFmanager.cxx
 * @brief Implementation of DSL::ASDFmanager: discovery of the HDUs of an ASDF tree.
 */
#include "ASDFmanager.h"

#include <algorithm>
#include <cctype>
#include <mutex>
#include <vector>

#include "ASDFheader.h"
#include "ASDFhandle.h"
#include "ASDFimage.h"
#include "ASDFsource.h"
#include "ASDFstate.h"
#include "ASDFtable.h"
#include "ASDFtree.h"

namespace DSL
{
    using namespace asdfdetail;

    namespace
    {
        const char* const kWho = "ASDFmanager"; //!< Name given in the exceptions thrown here

        /**
         * @brief Put a text in lower case, to compare HDU names.
         * @param text Text to convert.
         * @return The converted text.
         */
        std::string lower(std::string text)
        {
            std::transform(text.begin(), text.end(), text.begin(), [](unsigned char c){ return std::tolower(c); });
            return text;
        }

        /** @brief What the manager knows about one HDU; the data stay in the file until requested. */
        struct HDUEntry
        {
            /** @brief What an HDU holds. */
            enum Kind
            {
                empty,  //!< Only a header
                image,  //!< A numeric ndarray
                table   //!< A structured ndarray or a table value
            } kind = empty;                     //!< What this HDU holds
            std::string name;                   //!< EXTNAME or tree path, used to find the HDU
            FITShdu header;                     //!< Keywords of the HDU
            NdarrayPtr ndarray;                 //!< Image pixels, or the rows of a structured table
            std::vector<NdarrayPtr> arrays;     //!< One per column of a column-based table
            std::vector<ColumnSpec> columns;    //!< Layout of a table's columns
            uint64_t rows = 0;                  //!< Number of rows of a table
            bool structured = false;            //!< The table's columns all live in @c ndarray

            std::string dataPath;               //!< Path of the data in the tree
            std::string parentPath;             //!< Path of the node holding the data; empty at the root
            std::string keywordsPath;           //!< Path of the mapping whose scalars are keywords of the HDU
            bool hasKeywords = false;           //!< @c keywordsPath is set; false for data at the root
            std::string cardsPath;              //!< Astropy FITS layout: path of the card list
        };

        /**
         * @brief Tell whether a value is tagged as a table.
         * @param value Value of the ASDF tree.
         * @return True for a .../table/table-* or .../core/table-* tag.
         */
        bool isTableValue(asdf_value_t* value)
        {
            return hasTag(value, "table/table-") || hasTag(value, "core/table-");
        }

        /**
         * @brief Tell whether a value holds, at any depth, something that becomes an HDU.
         * @param value Value of the ASDF tree.
         * @return True if it is or contains data outside any GWCS object or transform.
         */
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

        /**
         * @brief Fill an HDU from its data, reading only metadata.
         * @details A table value or a structured ndarray makes a table. Any other ndarray with dimensions and
         * elements makes an image, whose type is only checked on loading; otherwise the HDU stays empty.
         * @param value Table or ndarray value.
         * @param entry HDU to fill; its header receives the structure keywords.
         * @throws FITSexception BAD_DIMEN if the columns of a table have different lengths, BAD_DATATYPE if the
         * ndarray can't be read, NOT_TABLE or BAD_TFORM if the table can't be described.
         */
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

    /** @brief Private state of DSL::ASDFmanager: the open file, its tree and its list of HDUs. */
    struct ASDFmanager::Impl
    {
        mutable std::mutex mtx;                 //!< Taken by every public method; before the file's own mutex
        std::shared_ptr<FileSource> source;     //!< The open file, shared with the ndarray handles; null once closed
        std::string path;                       //!< Path of the file, for messages
        ASDFnode tree;                          //!< The whole tree; its ndarray handles keep the file open
        std::vector<HDUEntry> hdus;             //!< All HDUs, the primary one first
        bool fitsLayout = false;                //!< The file holds the astropy FITS layout
        FITShdu rootKeywords;                   //!< Keywords from the scalars at the root of the tree

        /** @brief Release the ndarrays, which use the file, while nothing else does. */
        ~Impl()
        {
            if(source)
            {
                std::lock_guard<std::mutex> lock(source->mtx);
                hdus.clear();
            }
        }

        /**
         * @brief Add every HDU of a value in the astropy FITS layout.
         * @details An HDU is named after EXTNAME, otherwise PRIMARY for the first and HDUn for the others.
         * @param value    Sequence of {header, data} mappings; a non-sequence or non-mapping entry is skipped.
         * @param nodePath Path of the value in the tree.
         * @throws FITSexception As describeData().
         */
        void addFitsLayout(asdf_value_t* value, const std::string& nodePath)
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
                entry.parentPath = ASDFnode::ChildPath(nodePath, std::to_string(it->index));
                entry.cardsPath  = ASDFnode::ChildPath(entry.parentPath, "header");
                entry.dataPath   = ASDFnode::ChildPath(entry.parentPath, "data");
                if(const ASDFnode* cards = tree.FindPath(entry.cardsPath))
                    readCards(*cards, entry.header);

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

        /**
         * @brief Add one HDU for a table or an ndarray found in a tree of another layout.
         * @param value    Data of the HDU.
         * @param name     Tree path joined by dots, used as the name of the HDU and, unless it is PRIMARY, as EXTNAME.
         * @param nodePath Path of the value in the tree.
         * @param siblings Path of the mapping containing the data, whose scalars become keywords; null at the root.
         * @throws FITSexception As describeData().
         */
        void addData(asdf_value_t* value, const std::string& name, const std::string& nodePath, const std::string* siblings)
        {
            HDUEntry entry;
            entry.name = name;
            entry.dataPath = nodePath;
            entry.parentPath = ASDFnode::ParentPath(nodePath);
            entry.hasKeywords = siblings != nullptr;
            if(siblings)
            {
                entry.keywordsPath = *siblings;
                if(const ASDFnode* mapping = tree.FindPath(*siblings))
                    readTreeKeywords(*mapping, std::string(), *siblings, entry.header);
            }
            if(lower(name) != "primary")
                entry.header.ValueForKey("EXTNAME", name, fChar);
            describeData(value, entry);
            hdus.push_back(std::move(entry));
        }

        /**
         * @brief Find the HDUs below a value, recursively.
         * @details GWCS objects and transforms are skipped, as are asdf_library and history at the root.
         * @param value    Value to search.
         * @param name     Tree path of the value, its parts joined by dots; empty at the root.
         * @param nodePath Path of the value in the tree.
         * @param siblings Path of the mapping containing @p value, whose scalars become keywords if @p value is
         * data; null at the root and for its direct entries.
         * @throws FITSexception As describeData().
         */
        void walk(asdf_value_t* value, const std::string& name, const std::string& nodePath, const std::string* siblings)
        {
            if(isOpaqueValue(value))
                return;
            if(hasTag(value, "/fits/fits-"))
            {
                addFitsLayout(value, nodePath);
                return;
            }
            if(isDataValue(value))
            {
                addData(value, name, nodePath, siblings);
                return;
            }

            const std::string prefix = name.empty() ? std::string() : name + ".";
            asdf_mapping_t* mapping = nullptr;
            asdf_sequence_t* sequence = nullptr;
            if(asdf_value_as_mapping(value, &mapping) == ASDF_VALUE_OK)
            {
                asdf_mapping_iter_t* it = asdf_mapping_iter_init(mapping);
                while(asdf_mapping_iter_next(&it))
                {
                    const std::string key(it->key);
                    if(name.empty() && (key == "asdf_library" || key == "history"))
                        continue;
                    if(holdsData(it->value))
                        walk(it->value, prefix + key, ASDFnode::ChildPath(nodePath, key), name.empty() ? nullptr : &nodePath);
                }
            }
            else if(asdf_value_as_sequence(value, &sequence) == ASDF_VALUE_OK)
            {
                asdf_sequence_iter_t* it = asdf_sequence_iter_init(sequence);
                while(asdf_sequence_iter_next(&it))
                    if(holdsData(it->value))
                        walk(it->value, prefix + std::to_string(it->index), ASDFnode::ChildPath(nodePath, std::to_string(it->index)), siblings);
            }
        }

        /**
         * @brief Read the tree, list the HDUs of the file and put the primary one first.
         * @details Except in the astropy FITS layout, a tree without data gets an empty PRIMARY HDU and the
         * root keywords go to the first HDU.
         * @throws FITSexception FILE_NOT_OPENED if the file has no tree, or as describeData().
         */
        void discover()
        {
            ValuePtr root(asdf_get_value(source->file.get(), ""));
            if(!root)
                throw FITSexception(FILE_NOT_OPENED, kWho, "open", "no ASDF tree in "+path+": "+lastError(source->file.get()));

            tree = readTree(root.get(), source, std::string());
            walk(root.get(), std::string(), std::string(), nullptr);

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

            readRootKeywords(tree, rootKeywords);
            hdus.front().header.Import(rootKeywords);
        }

        /**
         * @brief Check that the file is open.
         * @param what Name of the calling method, for the message.
         * @throws FITSexception FILE_NOT_OPENED if the file is closed.
         */
        void requireOpen(const char* what) const
        {
            if(!source)
                throw FITSexception(FILE_NOT_OPENED, kWho, what, "the file is closed");
        }

        /**
         * @brief Find an HDU from its index.
         * @param index Index of the HDU, starting at 1.
         * @param what  Name of the calling method, for the message.
         * @return The HDU.
         * @throws FITSexception FILE_NOT_OPENED if the file is closed, BAD_HDU_NUM if there is no such HDU.
         */
        HDUEntry& at(int index, const char* what)
        {
            requireOpen(what);
            if(index < 1 || index > static_cast<int>(hdus.size()))
                throw FITSexception(BAD_HDU_NUM, kWho, what, "HDU "+std::to_string(index)+" doesn't exist in "+path+" ("+std::to_string(hdus.size())+" HDU)");
            return hdus[static_cast<size_t>(index - 1)];
        }

        /**
         * @brief Find the first HDU of a name.
         * @param name Name of the HDU, matched without regard to case.
         * @param what Name of the calling method, for the message.
         * @return The HDU.
         * @throws FITSexception FILE_NOT_OPENED if the file is closed, BAD_HDU_NUM if there is no such HDU.
         */
        HDUEntry& named(const std::string& name, const char* what)
        {
            requireOpen(what);
            auto it = std::find_if(hdus.begin(), hdus.end(), [&](const HDUEntry& e){ return lower(e.name) == lower(name); });
            if(it == hdus.end())
                throw FITSexception(BAD_HDU_NUM, kWho, what, "no HDU named '"+name+"' in "+path);
            return *it;
        }

        /**
         * @brief Read the image of an HDU, or part of it, with the tree of the file.
         * @param entry HDU to read.
         * @param tile  Part to read; the whole image by default.
         * @return The image, an ASDFimg\<T\>.
         * @throws FITSexception NOT_IMAGE if the HDU is a table or empty, or as readImage().
         */
        std::shared_ptr<FITScube> loadImage(HDUEntry& entry, const Tile& tile = {})
        {
            if(entry.kind == HDUEntry::table)
                throw FITSexception(NOT_IMAGE, kWho, "GetImage", "HDU '"+entry.name+"' of "+path+" is a table");
            if(entry.kind == HDUEntry::empty)
                throw FITSexception(NOT_IMAGE, kWho, "GetImage", "HDU '"+entry.name+"' of "+path+" holds no data");

            std::shared_ptr<FITScube> cube;
            {
                std::lock_guard<std::mutex> lock(source->mtx);
                cube = readImage(entry.ndarray.get(), entry.header, tile, true);
            }

            ImageState state;
            state.tree = tree.Clone();
            state.imagePath = entry.dataPath;
            state.companions = !fitsLayout && !entry.parentPath.empty();
            state.link.cardsPath = entry.cardsPath;
            state.link.keywordsPath = entry.keywordsPath;
            state.link.hasKeywords = entry.hasKeywords;
            state.link.rootKeywords = !fitsLayout && &entry == &hdus.front();
            state.link.sinkPath = entry.hasKeywords ? entry.keywordsPath : std::string();
            CubeAccess::attach(dynamic_cast<ASDFcube&>(*cube), std::move(state));
            return cube;
        }

        /**
         * @brief Read the table of an HDU.
         * @param entry HDU to read.
         * @return The table.
         * @throws FITSexception NOT_TABLE if the HDU isn't a table, or as readTable().
         */
        std::shared_ptr<FITStable> loadTable(HDUEntry& entry)
        {
            if(entry.kind != HDUEntry::table)
                throw FITSexception(NOT_TABLE, kWho, "GetTable", "HDU '"+entry.name+"' of "+path+" isn't a table");
            std::lock_guard<std::mutex> lock(source->mtx);
            return readTable(entry.columns, entry.rows, entry.structured, entry.header);
        }
    };

    ASDFmanager::ASDFmanager(const std::string& path): ASDFmanager(path, ASDFreadOptions()) {}

    ASDFmanager::ASDFmanager(const std::string& path, const ASDFreadOptions& options): impl(std::make_unique<Impl>())
    {
        impl->path = path;

        auto source = std::make_shared<FileSource>();
        source->path = path;
        source->tmpDir = options.tmpDir;

        // Fields left at zero keep libasdf's defaults
        asdf_config_t config{};
        config.decomp.max_memory_bytes = options.maxMemoryBytes;
        if(!source->tmpDir.empty())
            config.decomp.tmp_dir = source->tmpDir.c_str();

        source->file.reset(asdf_open_file_ex(path.c_str(), "r", &config));
        if(!source->file)
            throw FITSexception(FILE_NOT_OPENED, kWho, "ctor", "can't open "+path+": "+lastError(nullptr));
        if(asdf_error(source->file.get()))
            throw FITSexception(FILE_NOT_OPENED, kWho, "ctor", "can't open "+path+": "+lastError(source->file.get()));

        impl->source = std::move(source);
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
        return static_cast<bool>(impl->source);
    }

    void ASDFmanager::Close()
    {
        std::lock_guard<std::mutex> lk(impl->mtx);
        if(!impl->source)
            return;

        // The ndarrays depend on the file: they go first, while nothing else uses it
        {
            std::lock_guard<std::mutex> lock(impl->source->mtx);
            impl->hdus.clear();
        }
        // Handles handed out keep the file open; the manager's own go now
        impl->tree = ASDFnode();
        impl->source.reset();
    }

    ASDFnode ASDFmanager::GetTree()
    {
        std::lock_guard<std::mutex> lk(impl->mtx);
        impl->requireOpen("GetTree");
        return impl->tree.Clone();
    }

    std::shared_ptr<ASDFarray> ASDFmanager::GetArray(const std::string& path)
    {
        std::lock_guard<std::mutex> lk(impl->mtx);
        impl->requireOpen("GetArray");
        const ASDFnode* node = impl->tree.FindPath(path);
        if(!node)
            throw FITSexception(KEY_NO_EXIST, kWho, "GetArray", "nothing at '"+path+"' in "+impl->path);
        if(!node->IsArray())
            throw FITSexception(NOT_IMAGE, kWho, "GetArray", "'"+path+"' of "+impl->path+" isn't an ndarray");
        return node->Array()->Clone();
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

    ASDFhduPath ASDFmanager::GetPathAtIndex(const int& index)
    {
        std::lock_guard<std::mutex> lk(impl->mtx);
        const HDUEntry& entry = impl->at(index, "GetPathAtIndex");

        ASDFhduPath where;
        // Only the entries of the astropy FITS layout have a card list; keywordsPath is empty without keywords
        where.node = entry.cardsPath.empty() ? entry.keywordsPath : entry.parentPath;
        where.data = entry.dataPath;
        return where;
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

    std::shared_ptr<FITScube> ASDFmanager::GetImageWindow(const std::string& name, const std::vector<size_t>& origin, const std::vector<size_t>& size)
    {
        std::lock_guard<std::mutex> lk(impl->mtx);
        HDUEntry& entry = impl->named(name, "GetImageWindow");
        if(entry.kind != HDUEntry::image)
            return impl->loadImage(entry);     // throws the reason

        const std::vector<uint64_t> full = fitsAxes(entry.ndarray.get());
        const size_t rank = full.size();
        if(origin.size() > rank || size.size() > rank)
            throw FITSexception(BAD_DIMEN, kWho, "GetImageWindow", "HDU '"+name+"' has "+std::to_string(rank)+" axes");

        std::vector<size_t> start(origin);
        std::vector<size_t> length(size);
        start.resize(rank, 0);
        length.resize(rank, 0);
        Tile tile;
        tile.origin.assign(rank, 0);
        tile.size.assign(rank, 0);
        for(size_t k = 0; k < rank; k++)
        {
            if(start[k] >= full[k])
                throw FITSexception(BAD_DIMEN, kWho, "GetImageWindow", "position "+std::to_string(start[k])+" is past the end of axis "+std::to_string(k + 1)+" ("+std::to_string(full[k])+" pixels)");
            if(length[k] == 0 || start[k] + length[k] > full[k])
                length[k] = full[k] - start[k];
            tile.origin[rank - 1 - k] = start[k];
            tile.size  [rank - 1 - k] = length[k];
        }

        std::shared_ptr<FITScube> cube = impl->loadImage(entry, tile);
        CubeAccess::afterWindow(dynamic_cast<ASDFcube&>(*cube), std::vector<size_t>(full.begin(), full.end()), start, "GetImageWindow");
        return cube;
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
