//
//  ASDFheader.cxx
//  asdf2fits
//
//  Written by Claude (Anthropic's AI model) for William Gillard.
//  Copyright (c) 2026 William Gillard.
//
//  SPDX-License-Identifier: BSD-3-Clause
//  Licence text: see the LICENSE file at the root of this project.

/**
 * @file ASDFheader.cxx
 * @brief Implementation of the conversion between the ASDF tree and FITShdu keywords.
 */
#include "ASDFheader.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <regex>
#include <set>

namespace DSL::asdfdetail
{
    namespace
    {
        /**
         * @brief Put a text in upper case, to name keywords.
         * @param text Text to convert.
         * @return The converted text.
         */
        std::string upper(std::string text)
        {
            std::transform(text.begin(), text.end(), text.begin(), [](unsigned char c){ return std::toupper(c); });
            return text;
        }

        /**
         * @brief Put a text in lower case, to name tree keys.
         * @param text Text to convert.
         * @return The converted text.
         */
        std::string lower(std::string text)
        {
            std::transform(text.begin(), text.end(), text.begin(), [](unsigned char c){ return std::tolower(c); });
            return text;
        }

        /**
         * @brief Remove the spaces and tabs around a text.
         * @param text Text to trim.
         * @return The trimmed text; empty if nothing else is left.
         */
        std::string trim(const std::string& text)
        {
            const size_t first = text.find_first_not_of(" \t");
            if(first == std::string::npos)
                return std::string();
            return text.substr(first, text.find_last_not_of(" \t") - first + 1);
        }

        /**
         * @brief Add text to the single entry FITShdu keeps for a commentary keyword.
         * @param hdu  Header to update.
         * @param key  COMMENT or HISTORY.
         * @param text Text to add, after a space if the entry already exists.
         */
        void appendCommentary(FITShdu& hdu, const std::string& key, const std::string& text)
        {
            std::string all = hdu.Exists(key) ? hdu.GetEntry(key)->second.value() + " " + text : text;
            hdu.ValueForKey(key, all, fChar);
        }

        /**
         * @brief Set a keyword with the type of a scalar node.
         * @details A null node, a container or an ndarray sets nothing. Infinities and NaN are kept as text.
         * @param value   Scalar node.
         * @param key     Name of the keyword.
         * @param comment Comment of the keyword.
         * @param hdu     Header receiving the keyword.
         * @return True if a keyword was set.
         */
        bool setScalar(const ASDFnode& value, const std::string& key, const std::string& comment, FITShdu& hdu)
        {
            switch(value.GetKind())
            {
                case ASDFnode::Kind::boolean:
                    hdu.ValueForKey(key, value.Get<bool>(), comment);
                    return true;
                case ASDFnode::Kind::unsigned_integer:
                    hdu.ValueForKey(key, value.Get<uint64_t>(), comment);
                    return true;
                case ASDFnode::Kind::integer:
                    hdu.ValueForKey(key, value.Get<int64_t>(), comment);
                    return true;
                case ASDFnode::Kind::real:
                    // FITS has no form for infinities and NaN: they stay as text
                    if(std::isfinite(value.Get<double>()))
                        hdu.ValueForKey(key, value.Get<double>(), comment);
                    else
                        hdu.ValueForKey(key, value.Text(), fChar, comment);
                    return true;
                case ASDFnode::Kind::string:
                    hdu.ValueForKey(key, value.Get<std::string>(), fChar, comment);
                    return true;
                default:
                    return false;
            }
        }

        /**
         * @brief Get the key of a card.
         * @param card Card of an astropy FITS header.
         * @param key  Receives its key.
         * @return False for a card without a text key, such as a blank card.
         */
        bool cardKey(const ASDFnode& card, std::string& key)
        {
            if(!card.IsSequence() || card.Size() < 1 || !card[0].IsString())
                return false;
            key = card[0].Get<std::string>();
            return true;
        }

        /**
         * @brief Get the value of a keyword as a node of its keyword type.
         * @details A value that doesn't parse as its declared type is kept as text.
         * @param hdu Header holding the keyword.
         * @param key Name of the keyword, which exists.
         * @return The value.
         */
        ASDFnode keywordValue(const FITShdu& hdu, const std::string& key)
        {
            const FITSkeyword& keyword = hdu.GetEntry(key)->second;
            try
            {
                switch(keyword.type())
                {
                    case fBool:
                        return ASDFnode(hdu.GetBoolValueForKey(key));
                    case fShort: case fInt: case fLong: case fLongLong:
                        return ASDFnode(hdu.GetInt64ValueForKey(key));
                    case fUShort: case fUInt: case fULong: case fULongLong: case fByte:
                        return ASDFnode(hdu.GetUInt64ValueForKey(key));
                    case fFloat: case fDouble:
                        return ASDFnode(hdu.GetDoubleValueForKey(key));
                    default:
                        return ASDFnode(keyword.value());
                }
            }
            catch(const FITSexception&)
            {
                return ASDFnode(keyword.value());
            }
        }

        /**
         * @brief Make the [key, value, comment] card of a keyword.
         * @param hdu Header holding the keyword.
         * @param key Name of the keyword, which exists.
         * @return The card, without a comment if the keyword has none.
         */
        ASDFnode cardOf(const FITShdu& hdu, const std::string& key)
        {
            ASDFnode card = ASDFnode::Sequence();
            card.Append(key);
            card.Append(keywordValue(hdu, key));
            const std::string comment = trim(hdu.GetEntry(key)->second.comment());
            if(!comment.empty())
                card.Append(comment);
            return card;
        }

        /**
         * @brief Tell whether a keyword is the same in two headers.
         * @param a   First header.
         * @param b   Second header.
         * @param key Name of the keyword.
         * @return True if both have it with the same value, type and comment, or neither has it.
         */
        bool sameKeyword(const FITShdu& a, const FITShdu& b, const std::string& key)
        {
            auto ia = a.GetEntry(key);
            auto ib = b.GetEntry(key);
            const bool inA = ia != a.end();
            const bool inB = ib != b.end();
            if(inA != inB)
                return false;
            return !inA || (ia->second.value() == ib->second.value() && ia->second.type() == ib->second.type()
                            && ia->second.comment() == ib->second.comment());
        }

        /**
         * @brief Remove the cards of a keyword.
         * @param cards Card list.
         * @param key   Keyword; for COMMENT, blank-key cards holding text go too, as readCards() counts them as COMMENT.
         */
        void eraseCards(ASDFnode& cards, const std::string& key)
        {
            std::vector<ASDFnode>& items = cards.Items();
            items.erase(std::remove_if(items.begin(), items.end(), [&](const ASDFnode& card)
            {
                std::string k;
                if(!cardKey(card, k))
                    return false;
                return k == key || (key == "COMMENT" && k.empty() && card.Size() > 1 && !card[1].Text().empty());
            }), items.end());
        }

        /**
         * @brief Append a commentary keyword as cards of 72 characters of text.
         * @param cards Card list.
         * @param key   COMMENT or HISTORY.
         * @param text  Whole text of the keyword.
         */
        void appendCommentaryCards(ASDFnode& cards, const std::string& key, const std::string& text)
        {
            const size_t width = 72;
            for(size_t pos = 0; pos < text.size(); pos += width)
            {
                ASDFnode card = ASDFnode::Sequence();
                card.Append(key);
                card.Append(text.substr(pos, width));
                cards.Append(std::move(card));
            }
        }

        /**
         * @brief Copy into an astropy card list the changes made to a header.
         * @param cards    Card list to update.
         * @param hdu      The header now.
         * @param snapshot The header when it matched @p cards.
         */
        void syncCards(ASDFnode& cards, const FITShdu& hdu, const FITShdu& snapshot)
        {
            if(!cards.IsSequence())
                cards = ASDFnode::Sequence();

            // Removed keywords
            for(auto it = snapshot.begin(); it != snapshot.end(); ++it)
                if(!hdu.Exists(it->first))
                    eraseCards(cards, it->first);

            // Changed and new keywords
            for(auto it = hdu.begin(); it != hdu.end(); ++it)
            {
                const std::string& key = it->first;
                if(sameKeyword(hdu, snapshot, key))
                    continue;

                if(key == "COMMENT" || key == "HISTORY")
                {
                    eraseCards(cards, key);
                    appendCommentaryCards(cards, key, it->second.value());
                    continue;
                }

                bool replaced = false;
                for(ASDFnode& card : cards.Items())
                {
                    std::string k;
                    if(cardKey(card, k) && k == key)
                    {
                        card = cardOf(hdu, key);
                        replaced = true;
                        break;
                    }
                }
                if(!replaced)
                    cards.Append(cardOf(hdu, key));
            }
        }

        /**
         * @brief Find where a keyword with no place in the tree goes.
         * @details A dotted name A.B goes to the nested mappings a/b below the sink, if they can be made;
         * any other name goes to the sink as it is, so that reading the tree again gives the same keyword.
         * @param tree Tree to update.
         * @param sink Path of the sink mapping.
         * @param key  Name of the keyword.
         * @return The path to set.
         */
        std::string newKeywordPath(const ASDFnode& tree, const std::string& sink, const std::string& key)
        {
            if(key.find('.') == std::string::npos || key.front() == '.' || key.back() == '.')
                return ASDFnode::ChildPath(sink, key);

            std::string path = sink;
            size_t start = 0;
            while(true)
            {
                const size_t dot = key.find('.', start);
                const std::string part = lower(key.substr(start, dot == std::string::npos ? std::string::npos : dot - start));
                if(part.empty())
                    return ASDFnode::ChildPath(sink, key);
                path = ASDFnode::ChildPath(path, part);
                if(dot == std::string::npos)
                    break;

                // Every level but the last must be, or become, a mapping
                const ASDFnode* level = tree.FindPath(path);
                if(level && !level->IsMapping())
                    return ASDFnode::ChildPath(sink, key);
                start = dot + 1;
            }

            const ASDFnode* existing = tree.FindPath(path);
            return (existing && !existing->IsScalar()) ? ASDFnode::ChildPath(sink, key) : path;
        }
    }

#pragma region * Reading the tree
    bool isDataValue(asdf_value_t* value)
    {
        return hasTag(value, "core/ndarray-") || hasTag(value, "table/table-") || hasTag(value, "core/table-")
            || hasTag(value, "/fits/fits-");
    }

    bool isOpaqueValue(asdf_value_t* value)
    {
        return hasTag(value, "gwcs/") || hasTag(value, "/transform");
    }

    bool isDataNode(const ASDFnode& node)
    {
        return node.IsArray() || node.HasTag("core/ndarray-") || node.HasTag("table/table-") || node.HasTag("core/table-")
            || node.HasTag("/fits/fits-");
    }

    bool isOpaqueNode(const ASDFnode& node)
    {
        return node.HasTag("gwcs/") || node.HasTag("/transform");
    }

    bool holdsData(const ASDFnode& node)
    {
        if(isOpaqueNode(node))
            return false;
        if(isDataNode(node))
            return true;
        for(const ASDFnode& item : node.Items())
            if(holdsData(item))
                return true;
        return false;
    }

    void readCards(const ASDFnode& cards, FITShdu& hdu)
    {
        if(!cards.IsSequence())
            return;

        for(const ASDFnode& card : cards.Items())
        {
            std::string key;
            if(!cardKey(card, key))
                continue;   // astropy writes blank cards as []

            const ASDFnode* value = card.Size() > 1 ? &card[1] : nullptr;
            std::string comment;
            if(card.Size() > 2 && card[2].IsString())
                comment = card[2].Get<std::string>();

            if(key.empty() || key == "COMMENT" || key == "HISTORY")
            {
                // FITShdu keeps a single entry for commentary keywords, joining the cards as it
                // does when reading a FITS file. Text under a blank keyword is commentary too.
                if(value && value->IsScalar() && !value->Text().empty())
                    appendCommentary(hdu, key.empty() ? "COMMENT" : key, value->Text());
                continue;
            }

            if(value)
                setScalar(*value, key, comment, hdu);
            else
                hdu.ValueForKey(key, std::string(), fChar, comment);
        }
    }

    bool readScalar(const ASDFnode& value, const std::string& key, FITShdu& hdu)
    {
        return setScalar(value, key, std::string(), hdu);
    }

    void readTreeKeywords(const ASDFnode& mapping, const std::string& prefix, const std::string& path, FITShdu& hdu, KeyPaths* paths)
    {
        if(!mapping.IsMapping())
            return;

        for(size_t i = 0; i < mapping.Size(); i++)
        {
            const ASDFnode& child = mapping.Items()[i];
            const std::string& name = mapping.Keys()[i];
            if(isDataNode(child) || isOpaqueNode(child))
                continue;

            const std::string key = prefix.empty() ? upper(name) : prefix + "." + upper(name);
            const std::string childPath = ASDFnode::ChildPath(path, name);
            if(child.IsMapping())
                readTreeKeywords(child, key, childPath, hdu, paths);
            else if(readScalar(child, key, hdu) && paths)
                (*paths)[key] = childPath;
        }
    }

    void readRootKeywords(const ASDFnode& root, FITShdu& hdu, KeyPaths* paths)
    {
        if(!root.IsMapping())
            return;

        for(size_t i = 0; i < root.Size(); i++)
        {
            const ASDFnode& child = root.Items()[i];
            const std::string& name = root.Keys()[i];
            if(name == "asdf_library" || name == "history" || holdsData(child) || isOpaqueNode(child))
                continue;

            const std::string key = upper(name);
            const std::string childPath = ASDFnode::ChildPath(std::string(), name);
            if(child.IsMapping())
                readTreeKeywords(child, key, childPath, hdu, paths);
            else if(readScalar(child, key, hdu) && paths)
                (*paths)[key] = childPath;
        }
    }

#pragma endregion
#pragma region * Structure keywords
    void clearImageStructure(FITShdu& hdu)
    {
        std::vector<std::string> stale;
        static const std::regex naxis("NAXIS[0-9]*");
        for(auto it = hdu.begin(); it != hdu.end(); ++it)
            if(it->first == "BSCALE" || it->first == "BZERO" || std::regex_match(it->first, naxis))
                stale.push_back(it->first);
        for(const std::string& key : stale)
            hdu.DeleteKey(key);
    }

    void setImageStructure(FITShdu& hdu, asdf_scalar_datatype_t type, const std::vector<uint64_t>& naxis)
    {
        clearImageStructure(hdu);

        hdu.ValueForKey("BSCALE", 1.);
        switch(type)
        {
            case ASDF_DATATYPE_INT8:   hdu.ValueForKey("BITPIX", int32_t(8));   hdu.ValueForKey("BZERO", int32_t(-128));            break;
            case ASDF_DATATYPE_UINT8:
            case ASDF_DATATYPE_BOOL8:  hdu.ValueForKey("BITPIX", int32_t(8));   hdu.ValueForKey("BZERO", 0.);                       break;
            case ASDF_DATATYPE_INT16:  hdu.ValueForKey("BITPIX", int32_t(16));  hdu.ValueForKey("BZERO", 0.);                       break;
            case ASDF_DATATYPE_UINT16: hdu.ValueForKey("BITPIX", int32_t(16));  hdu.ValueForKey("BZERO", uint64_t(32768));          break;
            case ASDF_DATATYPE_INT32:  hdu.ValueForKey("BITPIX", int32_t(32));  hdu.ValueForKey("BZERO", 0.);                       break;
            case ASDF_DATATYPE_UINT32: hdu.ValueForKey("BITPIX", int32_t(32));  hdu.ValueForKey("BZERO", uint64_t(2147483648ULL));  break;
            case ASDF_DATATYPE_INT64:  hdu.ValueForKey("BITPIX", int32_t(64));  hdu.ValueForKey("BZERO", 0.);                       break;
            case ASDF_DATATYPE_UINT64: hdu.ValueForKey("BITPIX", int32_t(64));  hdu.ValueForKey("BZERO", uint64_t(1) << 63);        break;
            case ASDF_DATATYPE_FLOAT16:
            case ASDF_DATATYPE_FLOAT32:hdu.ValueForKey("BITPIX", int32_t(-32)); hdu.ValueForKey("BZERO", 0.);                       break;
            case ASDF_DATATYPE_FLOAT64:hdu.ValueForKey("BITPIX", int32_t(-64)); hdu.ValueForKey("BZERO", 0.);                       break;
            default:
                hdu.DeleteKey("BSCALE");
                hdu.ValueForKey("BITPIX", int32_t(8));
                break;
        }

        hdu.ValueForKey("NAXIS", static_cast<int32_t>(naxis.size()));
        for(size_t i = 0; i < naxis.size(); i++)
            hdu.ValueForKey("NAXIS" + std::to_string(i + 1), naxis[i]);
    }

    void setTableStructure(FITShdu& hdu, const std::vector<ColumnSpec>& columns, uint64_t rows)
    {
        static const std::regex stale("(NAXIS[0-9]*|BSCALE|BZERO|BLANK)");
        std::vector<std::string> remove;
        for(auto it = hdu.begin(); it != hdu.end(); ++it)
            if(std::regex_match(it->first, stale))
                remove.push_back(it->first);
        for(const std::string& key : remove)
            hdu.DeleteKey(key);

        uint64_t rowBytes = 0;
        for(const ColumnSpec& c : columns)
            rowBytes += c.elemBytes * c.repeat;

        hdu.ValueForKey("XTENSION", std::string("BINTABLE"), fChar, "binary table extension");
        hdu.ValueForKey("BITPIX", int32_t(8));
        hdu.ValueForKey("NAXIS", int32_t(2));
        hdu.ValueForKey("NAXIS1", rowBytes);
        hdu.ValueForKey("NAXIS2", rows);
        hdu.ValueForKey("PCOUNT", uint64_t(0));
        hdu.ValueForKey("GCOUNT", uint64_t(1));
        hdu.ValueForKey("TFIELDS", static_cast<int32_t>(columns.size()));

        for(size_t i = 0; i < columns.size(); i++)
        {
            const ColumnSpec& c = columns[i];
            const std::string n = std::to_string(i + 1);
            std::string form;
            switch(c.type)
            {
                case ASDF_DATATYPE_INT8:       form = "S"; break;
                case ASDF_DATATYPE_UINT8:      form = "B"; break;
                case ASDF_DATATYPE_INT16:      form = "I"; break;
                case ASDF_DATATYPE_UINT16:     form = "U"; break;
                case ASDF_DATATYPE_INT32:      form = "J"; break;
                case ASDF_DATATYPE_UINT32:     form = "V"; break;
                case ASDF_DATATYPE_INT64:      form = "K"; break;
                case ASDF_DATATYPE_UINT64:     form = "W"; break;
                case ASDF_DATATYPE_FLOAT32:    form = "E"; break;
                case ASDF_DATATYPE_FLOAT64:    form = "D"; break;
                case ASDF_DATATYPE_COMPLEX64:  form = "C"; break;
                case ASDF_DATATYPE_COMPLEX128: form = "M"; break;
                case ASDF_DATATYPE_BOOL8:      form = "L"; break;
                default:                       form = "A"; break;
            }

            const uint64_t count = (form == "A") ? c.elemBytes * c.repeat : c.repeat;
            hdu.ValueForKey("TTYPE" + n, c.name, fChar);
            hdu.ValueForKey("TFORM" + n, std::to_string(count) + form, fChar);
            if(!c.unit.empty())
                hdu.ValueForKey("TUNIT" + n, c.unit, fChar);
        }
    }

#pragma endregion
#pragma region * Writing cards with libasdf
    namespace
    {
        /**
         * @brief Append one [key, value, comment] card, with the value in its keyword type.
         * @details A value that doesn't parse as its declared type is written as text.
         * @param file    File the card is created in.
         * @param cards   List receiving the card.
         * @param key     Name of the keyword.
         * @param keyword Keyword to write; the comment is left out if it is blank.
         * @param hdu     Header holding the keyword, used to convert its value.
         */
        void appendCard(asdf_file_t* file, asdf_sequence_t* cards, const std::string& key, const FITSkeyword& keyword, const FITShdu& hdu)
        {
            asdf_sequence_t* card = asdf_sequence_create(file);
            asdf_sequence_set_style(card, ASDF_YAML_NODE_STYLE_FLOW);
            asdf_sequence_append_string0(card, key.c_str());

            try
            {
                switch(keyword.type())
                {
                    case fBool:
                        asdf_sequence_append_bool(card, hdu.GetBoolValueForKey(key));
                        break;
                    case fShort: case fInt: case fLong: case fLongLong:
                        asdf_sequence_append_int64(card, hdu.GetInt64ValueForKey(key));
                        break;
                    case fUShort: case fUInt: case fULong: case fULongLong: case fByte:
                        asdf_sequence_append_uint64(card, hdu.GetUInt64ValueForKey(key));
                        break;
                    case fFloat: case fDouble:
                        asdf_sequence_append_double(card, hdu.GetDoubleValueForKey(key));
                        break;
                    default:
                        asdf_sequence_append_string0(card, keyword.value().c_str());
                        break;
                }
            }
            catch(const FITSexception&)
            {
                // The value doesn't parse as its declared type: keep its text rather than lose it
                asdf_sequence_t* text = asdf_sequence_create(file);
                asdf_sequence_set_style(text, ASDF_YAML_NODE_STYLE_FLOW);
                asdf_sequence_append_string0(text, key.c_str());
                asdf_sequence_append_string0(text, keyword.value().c_str());
                asdf_sequence_destroy(card);
                card = text;
            }

            const std::string comment = trim(keyword.comment());
            if(!comment.empty())
                asdf_sequence_append_string0(card, comment.c_str());

            asdf_sequence_append_sequence(cards, card);
        }

        /**
         * @brief Append a commentary keyword as cards of 72 characters of text.
         * @param file  File the cards are created in.
         * @param cards List receiving the cards.
         * @param key   COMMENT or HISTORY.
         * @param text  Whole text of the keyword.
         */
        void appendCommentary(asdf_file_t* file, asdf_sequence_t* cards, const std::string& key, const std::string& text)
        {
            // One card holds 72 characters of text
            const size_t width = 72;
            for(size_t pos = 0; pos < text.size(); pos += width)
            {
                asdf_sequence_t* card = asdf_sequence_create(file);
                asdf_sequence_set_style(card, ASDF_YAML_NODE_STYLE_FLOW);
                asdf_sequence_append_string0(card, key.c_str());
                asdf_sequence_append_string0(card, text.substr(pos, width).c_str());
                asdf_sequence_append_sequence(cards, card);
            }
        }
    }

    asdf_sequence_t* writeCards(asdf_file_t* file, const FITShdu& hdu, const std::vector<std::string>& first)
    {
        asdf_sequence_t* cards = asdf_sequence_create(file);
        std::set<std::string> done(first.begin(), first.end());
        done.insert({"COMMENT", "HISTORY"});

        for(const std::string& key : first)
            if(hdu.Exists(key))
                appendCard(file, cards, key, hdu.GetEntry(key)->second, hdu);

        for(auto it = hdu.begin(); it != hdu.end(); ++it)
            if(done.find(it->first) == done.end())
                appendCard(file, cards, it->first, it->second, hdu);

        for(const char* key : {"COMMENT", "HISTORY"})
            if(hdu.Exists(key))
                appendCommentary(file, cards, key, hdu.GetEntry(key)->second.value());

        return cards;
    }

#pragma endregion
#pragma region * Header and tree in step
    bool isStructuralKey(const std::string& key, bool withExtname)
    {
        static const std::set<std::string> structural = {"SIMPLE", "XTENSION", "BITPIX", "NAXIS", "BSCALE", "BZERO", "BLANK",
                                                          "EXTEND", "PCOUNT", "GCOUNT", "CHECKSUM", "DATASUM"};
        static const std::regex naxis("NAXIS[0-9]+");
        return structural.count(key) || std::regex_match(key, naxis) || (withExtname && key == "EXTNAME");
    }

    void syncHeader(ASDFnode& tree, HeaderLink& link, const FITShdu& hdu, const FITShdu& snapshot)
    {
        if(!link.cardsPath.empty())
        {
            ASDFnode* cards = tree.FindPath(link.cardsPath);
            if(!cards)
                cards = &tree.SetPath(link.cardsPath, ASDFnode::Sequence());
            syncCards(*cards, hdu, snapshot);
            return;
        }

        // Removed keywords
        for(auto it = snapshot.begin(); it != snapshot.end(); ++it)
        {
            const std::string& key = it->first;
            if(hdu.Exists(key) || isStructuralKey(key, true))
                continue;
            auto mapped = link.keyPaths.find(key);
            if(mapped != link.keyPaths.end())
            {
                tree.ErasePath(mapped->second);
                link.keyPaths.erase(mapped);
            }
        }

        // Changed and new keywords
        for(auto it = hdu.begin(); it != hdu.end(); ++it)
        {
            const std::string& key = it->first;
            if(isStructuralKey(key, true) || sameKeyword(hdu, snapshot, key))
                continue;

            auto mapped = link.keyPaths.find(key);
            if(mapped != link.keyPaths.end())
            {
                try
                {
                    tree.SetPath(mapped->second, keywordValue(hdu, key));
                    continue;
                }
                catch(const FITSexception&)
                {
                    // The tree changed shape under the keyword: give it a new place below
                    link.keyPaths.erase(mapped);
                }
            }

            std::string path = newKeywordPath(tree, link.sinkPath, key);
            try
            {
                tree.SetPath(path, keywordValue(hdu, key));
            }
            catch(const FITSexception&)
            {
                path = ASDFnode::ChildPath(link.sinkPath, key);
                tree.SetPath(path, keywordValue(hdu, key));
            }
            link.keyPaths[key] = path;
        }
    }

    void removeCards(ASDFnode& cards, const std::string& key)
    {
        if(cards.IsSequence())
            eraseCards(cards, key);
    }

    void rebuildHeader(const ASDFnode& tree, HeaderLink& link, FITShdu& hdu)
    {
        const bool cards = !link.cardsPath.empty();

        FITShdu fresh;
        KeyPaths paths;
        if(cards)
        {
            if(const ASDFnode* list = tree.FindPath(link.cardsPath))
                readCards(*list, fresh);
        }
        else
        {
            if(link.hasKeywords)
                if(const ASDFnode* mapping = tree.FindPath(link.keywordsPath))
                    readTreeKeywords(*mapping, std::string(), link.keywordsPath, fresh, &paths);
            if(link.rootKeywords)
                readRootKeywords(tree, fresh, &paths);
        }

        std::vector<std::string> gone;
        for(auto it = hdu.begin(); it != hdu.end(); ++it)
            if(!isStructuralKey(it->first, !cards) && !fresh.Exists(it->first))
                gone.push_back(it->first);
        for(const std::string& key : gone)
            hdu.DeleteKey(key);

        for(auto it = fresh.begin(); it != fresh.end(); ++it)
            if(!isStructuralKey(it->first, !cards))
                hdu.ValueForKey(it->first, it->second.value(), it->second.type(), it->second.comment());

        if(!cards)
            link.keyPaths = std::move(paths);
    }
#pragma endregion
}
