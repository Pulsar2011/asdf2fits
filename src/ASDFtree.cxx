//
//  ASDFtree.cxx
//  asdf2fits
//
//  Written by Claude (Anthropic's AI model) for William Gillard.
//  Copyright (c) 2026 William Gillard.
//
//  SPDX-License-Identifier: BSD-3-Clause
//  Licence text: see the LICENSE file at the root of this project.

/**
 * @file ASDFtree.cxx
 * @brief Implementation of the conversion of libasdf values to ASDFnode trees, and of their writing as YAML.
 */
#include "ASDFtree.h"

#include <algorithm>
#include <cctype>
#include <set>

namespace DSL::asdfdetail
{
    namespace
    {
        /** @brief Start of the YAML core tags, which the kind of a node already tells. */
        const std::string kYamlCoreTags = "tag:yaml.org,2002:";
        /** @brief Start of the tags written with the "!" handle. */
        const std::string kAsdfTags = "tag:stsci.edu:asdf/";

        /**
         * @brief Get the tag of a value worth keeping.
         * @param value Value of the tree.
         * @return Its tag, or empty if it has none or a YAML core tag.
         */
        std::string tagOf(asdf_value_t* value)
        {
            const char* tag = asdf_value_tag(value);
            if(!tag || std::string(tag).rfind(kYamlCoreTags, 0) == 0)
                return std::string();
            return tag;
        }

        /**
         * @brief Convert a value; ndarrays become handles only if @p arrays is set.
         * @param value  The value.
         * @param file   The open file, for the handles.
         * @param path   Path of the value.
         * @param arrays Turn ndarrays into handles; false to convert their description like any value.
         * @return The node.
         */
        ASDFnode convert(asdf_value_t* value, const std::shared_ptr<FileSource>& file, const std::string& path, bool arrays)
        {
            const std::string tag = tagOf(value);

            if(arrays && tag.find("core/ndarray-") != std::string::npos)
            {
                ASDFnode desc = convert(value, file, path, false);
                return ASDFnode(ArrayAccess::make(describeArray(file, value, path, std::move(desc))));
            }

            ASDFnode node;
            asdf_mapping_t* mapping = nullptr;
            asdf_sequence_t* sequence = nullptr;
            if(asdf_value_as_mapping(value, &mapping) == ASDF_VALUE_OK)
            {
                node = ASDFnode::Mapping(tag);
                asdf_mapping_iter_t* it = asdf_mapping_iter_init(mapping);
                while(asdf_mapping_iter_next(&it))
                {
                    const std::string key(it->key);
                    node[key] = convert(it->value, file, ASDFnode::ChildPath(path, key), arrays);
                }
                return node;
            }
            if(asdf_value_as_sequence(value, &sequence) == ASDF_VALUE_OK)
            {
                node = ASDFnode::Sequence(tag);
                asdf_sequence_iter_t* it = asdf_sequence_iter_init(sequence);
                while(asdf_sequence_iter_next(&it))
                    node.Append(convert(it->value, file, ASDFnode::ChildPath(path, std::to_string(it->index)), arrays));
                return node;
            }

            bool     flag = false;
            int64_t  sint = 0;
            uint64_t uint = 0;
            double   real = 0.;
            const char* text = nullptr;

            // By exact type: asdf_value_as_bool() also accepts the integers 0 and 1
            switch(asdf_value_get_type(value))
            {
                case ASDF_VALUE_NULL:
                    break;
                case ASDF_VALUE_BOOL:
                    if(asdf_value_as_bool(value, &flag) == ASDF_VALUE_OK)
                        node = ASDFnode(flag);
                    break;
                case ASDF_VALUE_UINT64:
                    if(asdf_value_as_uint64(value, &uint) == ASDF_VALUE_OK)
                        node = ASDFnode(uint);
                    break;
                case ASDF_VALUE_INT8: case ASDF_VALUE_INT16: case ASDF_VALUE_INT32: case ASDF_VALUE_INT64:
                case ASDF_VALUE_UINT8: case ASDF_VALUE_UINT16: case ASDF_VALUE_UINT32:
                    if(asdf_value_as_int64(value, &sint) == ASDF_VALUE_OK)
                        node = ASDFnode(sint);
                    break;
                case ASDF_VALUE_FLOAT: case ASDF_VALUE_DOUBLE:
                    if(asdf_value_as_double(value, &real) == ASDF_VALUE_OK)
                        node = ASDFnode(real);
                    break;
                default:
                    if(asdf_value_as_scalar0(value, &text) == ASDF_VALUE_OK && text)
                        node = ASDFnode(std::string(text));
                    break;
            }

            // A value of a known type that didn't convert is kept as text
            if(node.IsNull() && asdf_value_get_type(value) != ASDF_VALUE_NULL
               && asdf_value_as_scalar0(value, &text) == ASDF_VALUE_OK && text)
                node = ASDFnode(std::string(text));

            node.SetTag(tag);
            return node;
        }

        /**
         * @brief Tell whether a string can be written without quotes and read back as the same string.
         * @param text The string.
         * @return True for names made of letters, digits, '_', '.', '-' and '/', starting with a letter or '_',
         * which aren't YAML words for null or booleans.
         */
        bool isPlain(const std::string& text)
        {
            if(text.empty())
                return false;
            const unsigned char first = static_cast<unsigned char>(text[0]);
            if(!std::isalpha(first) && first != '_')
                return false;
            for(unsigned char c : text)
                if(!std::isalnum(c) && c != '_' && c != '.' && c != '-' && c != '/')
                    return false;

            std::string word(text);
            std::transform(word.begin(), word.end(), word.begin(), [](unsigned char c){ return std::tolower(c); });
            static const std::set<std::string> reserved = {"null", "true", "false", "yes", "no", "on", "off", "y", "n"};
            return reserved.count(word) == 0;
        }

        /**
         * @brief Write a number in hexadecimal.
         * @param value  The number.
         * @param digits Number of digits.
         * @return The digits, in upper case.
         */
        std::string hex(uint32_t value, int digits)
        {
            static const char* const symbols = "0123456789ABCDEF";
            std::string text(static_cast<size_t>(digits), '0');
            for(int i = digits - 1; i >= 0; i--, value >>= 4)
                text[static_cast<size_t>(i)] = symbols[value & 0xF];
            return text;
        }

        /**
         * @brief Write a string in YAML double quotes.
         * @details Control characters, the YAML line breaks U+0085, U+2028 and U+2029, the byte order mark and
         * bytes that aren't valid UTF-8 are escaped, so that the document stays valid.
         * @param text The string.
         * @return The quoted string.
         */
        std::string doubleQuoted(const std::string& text)
        {
            std::string out = "\"";
            for(size_t i = 0; i < text.size();)
            {
                const unsigned char c = static_cast<unsigned char>(text[i]);
                if(c < 0x80)
                {
                    switch(c)
                    {
                        case '"':  out += "\\\""; break;
                        case '\\': out += "\\\\"; break;
                        case '\n': out += "\\n";  break;
                        case '\t': out += "\\t";  break;
                        case '\r': out += "\\r";  break;
                        default:
                            if(c < 0x20 || c == 0x7F)
                                out += "\\x" + hex(c, 2);
                            else
                                out += static_cast<char>(c);
                    }
                    i++;
                    continue;
                }

                // A UTF-8 sequence: its length from the first byte, then continuation bytes
                const size_t length = (c >= 0xC2 && c <= 0xDF) ? 2 : (c >= 0xE0 && c <= 0xEF) ? 3 : (c >= 0xF0 && c <= 0xF4) ? 4 : 0;
                bool valid = length > 0 && i + length <= text.size();
                uint32_t code = (length == 2) ? (c & 0x1F) : (length == 3) ? (c & 0x0F) : (c & 0x07);
                for(size_t k = 1; valid && k < length; k++)
                {
                    const unsigned char next = static_cast<unsigned char>(text[i + k]);
                    valid = (next & 0xC0) == 0x80;
                    code = (code << 6) | (next & 0x3F);
                }
                if(valid && length == 3 && (code < 0x800 || (code >= 0xD800 && code <= 0xDFFF)))
                    valid = false;
                if(valid && length == 4 && (code < 0x10000 || code > 0x10FFFF))
                    valid = false;

                if(!valid)
                {
                    out += "\\x" + hex(c, 2);
                    i++;
                    continue;
                }
                if(code < 0xA0 || code == 0x2028 || code == 0x2029 || code == 0xFEFF)
                    out += "\\u" + hex(code, 4);
                else
                    out.append(text, i, length);
                i += length;
            }
            return out + "\"";
        }

        /**
         * @brief Write a string as a scalar or a key.
         * @param text The string.
         * @return The text as it is, or quoted.
         */
        std::string stringText(const std::string& text)
        {
            return isPlain(text) ? text : doubleQuoted(text);
        }

        /**
         * @brief Write a tag.
         * @param tag Tag in full; may be empty.
         * @return "!short" for an ASDF tag, "!<tag>" for another one, followed by a space; empty for no tag.
         */
        std::string tagText(const std::string& tag)
        {
            if(tag.empty())
                return std::string();
            if(tag.rfind(kAsdfTags, 0) == 0 && tag.size() > kAsdfTags.size())
                return "!" + tag.substr(kAsdfTags.size()) + " ";
            return "!<" + tag + "> ";
        }

        /** @brief Writes a tree, line by line. */
        class Emitter
        {
        public:
            /** @brief Prepare to write. @param arrays Gives the node written in place of each ndarray. */
            explicit Emitter(const ArrayWriter& arrays): farrays(arrays) {}

            /**
             * @brief Write a node.
             * @param node   Node to write; an ndarray is replaced first.
             * @param path   Its path.
             * @param head   What comes before it on its line: "key:", "-" or "---".
             * @param indent Indentation of what is below it.
             */
            void write(const ASDFnode& node, const std::string& path, const std::string& head, int indent)
            {
                if(node.IsArray())
                {
                    write(farrays(path, node), path, head, indent);
                    return;
                }

                const std::string tag = tagText(node.Tag());
                if(node.IsScalar() || node.Size() == 0)
                {
                    fout += head + " " + tag + inlineText(node) + "\n";
                    return;
                }

                std::string line = head;
                if(!tag.empty())
                    line += " " + tag.substr(0, tag.size() - 1);
                fout += line + "\n";

                const std::string pad(static_cast<size_t>(indent), ' ');
                for(size_t i = 0; i < node.Size(); i++)
                {
                    const bool mapping = node.IsMapping();
                    const std::string part = mapping ? node.Keys()[i] : std::to_string(i);
                    const std::string itemHead = pad + (mapping ? stringText(part) + ":" : std::string("-"));
                    write(node.Items()[i], ASDFnode::ChildPath(path, part), itemHead, indent + 2);
                }
            }

            /** @brief Get what was written. @return The text. */
            std::string& text() { return fout; }

        private:
            /**
             * @brief Write a scalar, or an empty container, on one line.
             * @param node The node.
             * @return Its text, without its tag.
             */
            static std::string inlineText(const ASDFnode& node)
            {
                switch(node.GetKind())
                {
                    case ASDFnode::Kind::null:     return "null";
                    case ASDFnode::Kind::string:   return stringText(node.Get<std::string>());
                    case ASDFnode::Kind::sequence: return "[]";
                    case ASDFnode::Kind::mapping:  return "{}";
                    default:                       return node.Text();
                }
            }

            const ArrayWriter& farrays; //!< Gives the node written in place of each ndarray
            std::string fout;           //!< The document written so far
        };
    }

    ASDFnode readTree(asdf_value_t* value, const std::shared_ptr<FileSource>& file, const std::string& path)
    {
        return convert(value, file, path, true);
    }

    std::string writeDocument(const ASDFnode& root, const ArrayWriter& arrays)
    {
        Emitter emitter(arrays);
        emitter.text() = "#ASDF 1.0.0\n#ASDF_STANDARD 1.6.0\n%YAML 1.1\n%TAG ! tag:stsci.edu:asdf/\n";

        if(root.Tag().empty() && root.IsMapping())
        {
            ASDFnode tagged = root;
            tagged.SetTag(kRootTag);
            emitter.write(tagged, std::string(), "---", 0);
        }
        else
            emitter.write(root, std::string(), "---", 0);

        emitter.text() += "...\n";
        return std::move(emitter.text());
    }
}
