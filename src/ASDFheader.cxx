//
//  ASDFheader.cxx
//  asdf2fits
//
//  Written by Claude (Anthropic's AI model) for William Gillard.
//  Copyright (c) 2026 William Gillard.
//
//  SPDX-License-Identifier: BSD-3-Clause
//  Licence text: see the LICENSE file at the root of this project.

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
        std::string upper(std::string text)
        {
            std::transform(text.begin(), text.end(), text.begin(), [](unsigned char c){ return std::toupper(c); });
            return text;
        }

        std::string trim(const std::string& text)
        {
            const size_t first = text.find_first_not_of(" \t");
            if(first == std::string::npos)
                return std::string();
            return text.substr(first, text.find_last_not_of(" \t") - first + 1);
        }

        void appendCommentary(FITShdu& hdu, const std::string& key, const std::string& text)
        {
            std::string all = hdu.Exists(key) ? hdu.GetEntry(key)->second.value() + " " + text : text;
            hdu.ValueForKey(key, all, fChar);
        }

        void setScalar(asdf_value_t* value, const std::string& key, const std::string& comment, FITShdu& hdu)
        {
            bool     flag = false;
            int64_t  sint = 0;
            uint64_t uint = 0;
            double   real = 0.;
            const char* text = nullptr;

            // By exact type: asdf_value_as_bool() also accepts the integers 0 and 1
            switch(asdf_value_get_type(value))
            {
                case ASDF_VALUE_NULL:
                    return;
                case ASDF_VALUE_BOOL:
                    if(asdf_value_as_bool(value, &flag) == ASDF_VALUE_OK)
                        hdu.ValueForKey(key, flag, comment);
                    return;
                case ASDF_VALUE_UINT64:
                    if(asdf_value_as_uint64(value, &uint) == ASDF_VALUE_OK)
                    {
                        hdu.ValueForKey(key, uint, comment);
                        return;
                    }
                    break;
                case ASDF_VALUE_INT8: case ASDF_VALUE_INT16: case ASDF_VALUE_INT32: case ASDF_VALUE_INT64:
                case ASDF_VALUE_UINT8: case ASDF_VALUE_UINT16: case ASDF_VALUE_UINT32:
                    if(asdf_value_as_int64(value, &sint) == ASDF_VALUE_OK)
                    {
                        hdu.ValueForKey(key, sint, comment);
                        return;
                    }
                    break;
                case ASDF_VALUE_FLOAT: case ASDF_VALUE_DOUBLE:
                    // FITS has no form for infinities and NaN: they stay as text below
                    if(asdf_value_as_double(value, &real) == ASDF_VALUE_OK && std::isfinite(real))
                    {
                        hdu.ValueForKey(key, real, comment);
                        return;
                    }
                    break;
                default:
                    break;
            }

            if(asdf_value_as_scalar0(value, &text) == ASDF_VALUE_OK && text)
                hdu.ValueForKey(key, std::string(text), fChar, comment);
        }
    }

    bool isDataValue(asdf_value_t* value)
    {
        return hasTag(value, "core/ndarray-") || hasTag(value, "table/table-") || hasTag(value, "core/table-")
            || hasTag(value, "/fits/fits-");
    }

    bool isOpaqueValue(asdf_value_t* value)
    {
        return hasTag(value, "gwcs/") || hasTag(value, "/transform");
    }

    void readCards(asdf_value_t* cards, FITShdu& hdu)
    {
        asdf_sequence_t* list = nullptr;
        if(asdf_value_as_sequence(cards, &list) != ASDF_VALUE_OK)
            return;

        asdf_sequence_iter_t* it = asdf_sequence_iter_init(list);
        while(asdf_sequence_iter_next(&it))
        {
            asdf_sequence_t* card = nullptr;
            if(asdf_value_as_sequence(it->value, &card) != ASDF_VALUE_OK)
                continue;

            const int size = asdf_sequence_size(card);
            if(size < 1)
                continue;   // astropy writes blank cards as []

            ValuePtr keyValue(asdf_sequence_get(card, 0));
            const char* keyText = nullptr;
            if(!keyValue || asdf_value_as_string0(keyValue.get(), &keyText) != ASDF_VALUE_OK || !keyText)
                continue;
            const std::string key(keyText);

            ValuePtr value(size > 1 ? asdf_sequence_get(card, 1) : nullptr);
            ValuePtr note (size > 2 ? asdf_sequence_get(card, 2) : nullptr);

            const char* noteText = nullptr;
            std::string comment;
            if(note && asdf_value_as_string0(note.get(), &noteText) == ASDF_VALUE_OK && noteText)
                comment = noteText;

            const char* text = nullptr;
            if(key.empty() || key == "COMMENT" || key == "HISTORY")
            {
                // FITShdu keeps a single entry for commentary keywords, joining the cards as it
                // does when reading a FITS file. Text under a blank keyword is commentary too.
                if(value && asdf_value_as_scalar0(value.get(), &text) == ASDF_VALUE_OK && text && *text)
                    appendCommentary(hdu, key.empty() ? "COMMENT" : key, text);
                continue;
            }

            if(value)
                setScalar(value.get(), key, comment, hdu);
            else
                hdu.ValueForKey(key, std::string(), fChar, comment);
        }
    }

    void readScalar(asdf_value_t* value, const std::string& key, FITShdu& hdu)
    {
        if(!asdf_value_is_container(value))
            setScalar(value, key, std::string(), hdu);
    }

    void readTreeKeywords(asdf_value_t* value, const std::string& prefix, FITShdu& hdu)
    {
        asdf_mapping_t* mapping = nullptr;
        if(asdf_value_as_mapping(value, &mapping) != ASDF_VALUE_OK)
            return;

        asdf_mapping_iter_t* it = asdf_mapping_iter_init(mapping);
        while(asdf_mapping_iter_next(&it))
        {
            if(isDataValue(it->value) || isOpaqueValue(it->value))
                continue;

            const std::string key = prefix.empty() ? upper(it->key) : prefix + "." + upper(it->key);
            if(asdf_value_is_mapping(it->value))
                readTreeKeywords(it->value, key, hdu);
            else
                readScalar(it->value, key, hdu);
        }
    }

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

    namespace
    {
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
}
