//
//  ASDFnode.cxx
//  asdf2fits
//
//  Written by Claude (Anthropic's AI model) for William Gillard.
//  Copyright (c) 2026 William Gillard.
//
//  SPDX-License-Identifier: BSD-3-Clause
//  Licence text: see the LICENSE file at the root of this project.

/**
 * @file ASDFnode.cxx
 * @brief Implementation of DSL::ASDFnode.
 */
#include "ASDFnode.h"

#include <charconv>
#include <cmath>
#include <ostream>

#include <fitsio.h>

#include <DSTfits/FITSexception.h>

#include "ASDFarray.h"

namespace DSL
{
    namespace
    {
        const char* const kWho = "ASDFnode";    //!< Name given in the exceptions thrown here

        /**
         * @brief Name a kind of node, for messages.
         * @param kind Kind to name.
         * @return Its name.
         */
        const char* kindName(ASDFnode::Kind kind)
        {
            switch(kind)
            {
                case ASDFnode::Kind::null:             return "null";
                case ASDFnode::Kind::boolean:          return "a boolean";
                case ASDFnode::Kind::integer:          return "an integer";
                case ASDFnode::Kind::unsigned_integer: return "an unsigned integer";
                case ASDFnode::Kind::real:             return "a real";
                case ASDFnode::Kind::string:           return "a string";
                case ASDFnode::Kind::sequence:         return "a sequence";
                case ASDFnode::Kind::mapping:          return "a mapping";
                case ASDFnode::Kind::array:            return "an ndarray";
            }
            return "unknown";
        }

        /**
         * @brief Read a sequence index from a path part.
         * @param part   Path part.
         * @param result Receives the index.
         * @return False if the part isn't a plain decimal number.
         */
        bool parseIndex(const std::string& part, size_t& result)
        {
            if(part.empty())
                return false;
            const auto [end, ec] = std::from_chars(part.data(), part.data() + part.size(), result);
            return ec == std::errc() && end == part.data() + part.size();
        }
    }

#pragma region * ctor
    ASDFnode::ASDFnode(const char* text)
    {
        if(text)
        {
            fkind = Kind::string;
            ftext = text;
        }
    }

    ASDFnode::ASDFnode(std::shared_ptr<ASDFarray> array)
    {
        if(array)
        {
            fkind  = Kind::array;
            farray = std::move(array);
        }
    }

    ASDFnode ASDFnode::Mapping(std::string tag)
    {
        ASDFnode node;
        node.fkind = Kind::mapping;
        node.ftag  = std::move(tag);
        return node;
    }

    ASDFnode ASDFnode::Sequence(std::string tag)
    {
        ASDFnode node;
        node.fkind = Kind::sequence;
        node.ftag  = std::move(tag);
        return node;
    }

#pragma endregion
#pragma region * Scalars
    void ASDFnode::typeError(const char* wanted) const
    {
        throw FITSexception(BAD_DATATYPE, kWho, "Get", std::string("the node is ")+kindName(fkind)+", not "+wanted);
    }

    void ASDFnode::rangeError(const char* wanted) const
    {
        throw FITSexception(NUM_OVERFLOW, kWho, "Get", "the value "+Text()+" doesn't fit in "+wanted);
    }

    std::string ASDFnode::Text() const
    {
        switch(fkind)
        {
            case Kind::boolean:          return fbool ? "true" : "false";
            case Kind::integer:          return std::to_string(fint);
            case Kind::unsigned_integer: return std::to_string(fuint);
            case Kind::string:           return ftext;
            case Kind::real:
            {
                if(std::isnan(freal))
                    return ".nan";
                if(std::isinf(freal))
                    return freal > 0 ? ".inf" : "-.inf";

                // Shortest text that reads back as the same double, with a '.' so that YAML sees a real
                char buffer[64];
                const auto result = std::to_chars(buffer, buffer + sizeof(buffer), freal);
                std::string text(buffer, result.ptr);
                const size_t e = text.find_first_of("eE");
                std::string mantissa = text.substr(0, e);
                std::string exponent = (e == std::string::npos) ? std::string() : text.substr(e);
                if(mantissa.find('.') == std::string::npos)
                    mantissa += ".0";
                if(exponent.size() > 1 && exponent[1] != '+' && exponent[1] != '-')
                    exponent.insert(1, "+");
                return mantissa + exponent;
            }
            default:
                return std::string();
        }
    }

#pragma endregion
#pragma region * Containers
    const ASDFnode* ASDFnode::Find(const std::string& key) const
    {
        if(fkind != Kind::mapping)
            return nullptr;
        for(size_t i = 0; i < fkeys.size(); i++)
            if(fkeys[i] == key)
                return &fitems[i];
        return nullptr;
    }

    ASDFnode* ASDFnode::Find(const std::string& key)
    {
        return const_cast<ASDFnode*>(static_cast<const ASDFnode*>(this)->Find(key));
    }

    ASDFnode& ASDFnode::operator[](const std::string& key)
    {
        if(fkind == Kind::null)
            fkind = Kind::mapping;
        if(fkind != Kind::mapping)
            throw FITSexception(BAD_DATATYPE, kWho, "operator[]", std::string("can't look up '")+key+"' in "+kindName(fkind));

        if(ASDFnode* found = Find(key))
            return *found;
        fkeys.push_back(key);
        fitems.emplace_back();
        return fitems.back();
    }

    const ASDFnode& ASDFnode::operator[](const std::string& key) const
    {
        const ASDFnode* found = Find(key);
        if(!found)
            throw FITSexception(KEY_NO_EXIST, kWho, "operator[]", "no key '"+key+"' in "+kindName(fkind));
        return *found;
    }

    ASDFnode& ASDFnode::operator[](size_t index)
    {
        return const_cast<ASDFnode&>(static_cast<const ASDFnode&>(*this)[index]);
    }

    const ASDFnode& ASDFnode::operator[](size_t index) const
    {
        if(index >= fitems.size())
            throw FITSexception(BAD_DIMEN, kWho, "operator[]", "no item "+std::to_string(index)+" in "+kindName(fkind)+" of "+std::to_string(fitems.size())+" items");
        return fitems[index];
    }

    ASDFnode& ASDFnode::Append(ASDFnode value)
    {
        if(fkind == Kind::null)
            fkind = Kind::sequence;
        if(fkind != Kind::sequence)
            throw FITSexception(BAD_DATATYPE, kWho, "Append", std::string("can't append to ")+kindName(fkind));
        fitems.push_back(std::move(value));
        return fitems.back();
    }

    bool ASDFnode::Erase(const std::string& key)
    {
        if(fkind != Kind::mapping)
            return false;
        for(size_t i = 0; i < fkeys.size(); i++)
        {
            if(fkeys[i] == key)
            {
                fkeys.erase(fkeys.begin() + static_cast<std::ptrdiff_t>(i));
                fitems.erase(fitems.begin() + static_cast<std::ptrdiff_t>(i));
                return true;
            }
        }
        return false;
    }

    bool ASDFnode::EraseAt(size_t index)
    {
        if(index >= fitems.size())
            return false;
        if(fkind == Kind::mapping)
            fkeys.erase(fkeys.begin() + static_cast<std::ptrdiff_t>(index));
        fitems.erase(fitems.begin() + static_cast<std::ptrdiff_t>(index));
        return true;
    }

#pragma endregion
#pragma region * Paths
    std::string ASDFnode::JoinPath(const std::vector<std::string>& parts)
    {
        std::string path;
        for(const std::string& part : parts)
            path = ChildPath(path, part);
        return path;
    }

    std::string ASDFnode::ChildPath(const std::string& path, const std::string& part)
    {
        std::string escaped;
        for(char c : part)
        {
            if(c == '~')
                escaped += "~0";
            else if(c == '/')
                escaped += "~1";
            else
                escaped += c;
        }
        return path.empty() ? escaped : path + "/" + escaped;
    }

    std::vector<std::string> ASDFnode::SplitPath(const std::string& path)
    {
        std::vector<std::string> parts;
        if(path.empty())
            return parts;

        std::string part;
        for(size_t i = 0; i <= path.size(); i++)
        {
            if(i == path.size() || path[i] == '/')
            {
                parts.push_back(part);
                part.clear();
            }
            else if(path[i] == '~' && i + 1 < path.size() && (path[i+1] == '0' || path[i+1] == '1'))
            {
                part += (path[i+1] == '0') ? '~' : '/';
                i++;
            }
            else
                part += path[i];
        }
        return parts;
    }

    std::string ASDFnode::ParentPath(const std::string& path)
    {
        std::vector<std::string> parts = SplitPath(path);
        if(!parts.empty())
            parts.pop_back();
        return JoinPath(parts);
    }

    const ASDFnode* ASDFnode::FindPath(const std::string& path) const
    {
        const ASDFnode* node = this;
        for(const std::string& part : SplitPath(path))
        {
            if(node->fkind == Kind::mapping)
                node = node->Find(part);
            else if(node->fkind == Kind::sequence)
            {
                size_t index = 0;
                node = (parseIndex(part, index) && index < node->fitems.size()) ? &node->fitems[index] : nullptr;
            }
            else
                node = nullptr;

            if(!node)
                return nullptr;
        }
        return node;
    }

    ASDFnode* ASDFnode::FindPath(const std::string& path)
    {
        return const_cast<ASDFnode*>(static_cast<const ASDFnode*>(this)->FindPath(path));
    }

    const ASDFnode& ASDFnode::At(const std::string& path) const
    {
        const ASDFnode* node = FindPath(path);
        if(!node)
            throw FITSexception(KEY_NO_EXIST, kWho, "At", "nothing at '"+path+"'");
        return *node;
    }

    ASDFnode& ASDFnode::At(const std::string& path)
    {
        return const_cast<ASDFnode&>(static_cast<const ASDFnode&>(*this).At(path));
    }

    ASDFnode& ASDFnode::SetPath(const std::string& path, ASDFnode value)
    {
        const std::vector<std::string> parts = SplitPath(path);
        if(parts.empty())
            throw FITSexception(KEY_NO_EXIST, kWho, "SetPath", "the path is empty");

        ASDFnode* node = this;
        for(size_t i = 0; i < parts.size(); i++)
        {
            const std::string& part = parts[i];
            ASDFnode* next = nullptr;
            if(node->fkind == Kind::sequence)
            {
                size_t index = 0;
                if(!parseIndex(part, index) || index >= node->fitems.size())
                    throw FITSexception(KEY_NO_EXIST, kWho, "SetPath", "no item '"+part+"' in the sequence of '"+path+"'");
                next = &node->fitems[index];
            }
            else if(node->fkind == Kind::mapping || node->fkind == Kind::null)
                next = &(*node)[part];
            else
                throw FITSexception(BAD_DATATYPE, kWho, "SetPath", "'"+path+"' goes through "+kindName(node->fkind));

            // Intermediate nodes become mappings when they are missing
            if(i + 1 < parts.size() && next->fkind == Kind::null)
                next->fkind = Kind::mapping;
            node = next;
        }

        // A time or a unit replaced by another scalar stays tagged as such
        if(value.ftag.empty() && value.IsScalar() && !value.IsNull() && node->IsScalar() && !node->IsNull())
            value.ftag = node->ftag;
        *node = std::move(value);
        return *node;
    }

    bool ASDFnode::ErasePath(const std::string& path)
    {
        std::vector<std::string> parts = SplitPath(path);
        if(parts.empty())
            return false;

        const std::string last = parts.back();
        parts.pop_back();
        ASDFnode* parent = FindPath(JoinPath(parts));
        if(!parent)
            return false;
        if(parent->fkind == Kind::mapping)
            return parent->Erase(last);

        size_t index = 0;
        return parent->fkind == Kind::sequence && parseIndex(last, index) && parent->EraseAt(index);
    }

#pragma endregion
#pragma region * Tree
    void ASDFnode::Visit(const std::function<bool(const std::string&, ASDFnode&)>& visitor)
    {
        std::function<void(const std::string&, ASDFnode&)> walk = [&](const std::string& path, ASDFnode& node)
        {
            if(!visitor(path, node))
                return;
            for(size_t i = 0; i < node.fitems.size(); i++)
                walk(ChildPath(path, node.fkind == Kind::mapping ? node.fkeys[i] : std::to_string(i)), node.fitems[i]);
        };
        walk(std::string(), *this);
    }

    void ASDFnode::Visit(const std::function<bool(const std::string&, const ASDFnode&)>& visitor) const
    {
        std::function<void(const std::string&, const ASDFnode&)> walk = [&](const std::string& path, const ASDFnode& node)
        {
            if(!visitor(path, node))
                return;
            for(size_t i = 0; i < node.fitems.size(); i++)
                walk(ChildPath(path, node.fkind == Kind::mapping ? node.fkeys[i] : std::to_string(i)), node.fitems[i]);
        };
        walk(std::string(), *this);
    }

    ASDFnode ASDFnode::Clone() const
    {
        ASDFnode copy(*this);
        copy.Visit([](const std::string&, ASDFnode& node)
        {
            if(node.farray)
                node.farray = node.farray->Clone();
            return true;
        });
        return copy;
    }

    bool ASDFnode::operator==(const ASDFnode& other) const
    {
        if(ftag != other.ftag)
            return false;

        // libasdf reads small unsigned values as signed ones: compare integers by value
        if(IsInteger() && other.IsInteger())
        {
            if(fkind == other.fkind)
                return fkind == Kind::integer ? fint == other.fint : fuint == other.fuint;
            const int64_t  s = (fkind == Kind::integer) ? fint : other.fint;
            const uint64_t u = (fkind == Kind::integer) ? other.fuint : fuint;
            return std::cmp_equal(s, u);
        }

        if(fkind != other.fkind)
            return false;
        switch(fkind)
        {
            case Kind::null:             return true;
            case Kind::boolean:          return fbool == other.fbool;
            case Kind::integer:          return fint  == other.fint;
            case Kind::unsigned_integer: return fuint == other.fuint;
            case Kind::real:             return freal == other.freal || (std::isnan(freal) && std::isnan(other.freal));
            case Kind::string:           return ftext == other.ftext;
            case Kind::array:            return farray == other.farray;
            case Kind::sequence:
            case Kind::mapping:          return fkeys == other.fkeys && fitems == other.fitems;
        }
        return false;
    }

    void ASDFnode::Dump(std::ostream& out, int indent) const
    {
        const std::string pad(static_cast<size_t>(indent), ' ');
        const std::string tag = ftag.empty() ? std::string() : "!<" + ftag + "> ";

        auto inlineText = [](const ASDFnode& node) -> std::string
        {
            const std::string t = node.ftag.empty() ? std::string() : "!<" + node.ftag + "> ";
            switch(node.fkind)
            {
                case Kind::null:     return t + "null";
                case Kind::string:   return t + "\"" + node.ftext + "\"";
                case Kind::array:    return t + "<ndarray " + node.farray->DataType() + ">";
                case Kind::sequence: return node.fitems.empty() ? t + "[]" : std::string();
                case Kind::mapping:  return node.fitems.empty() ? t + "{}" : std::string();
                default:             return t + node.Text();
            }
        };

        if(fitems.empty())
        {
            out << pad << inlineText(*this) << "\n";
            return;
        }

        if(!tag.empty())
            out << pad << tag << "\n";
        for(size_t i = 0; i < fitems.size(); i++)
        {
            const ASDFnode& item = fitems[i];
            const std::string head = pad + (fkind == Kind::mapping ? fkeys[i] + ":" : std::string("-"));
            const std::string text = inlineText(item);
            if(!text.empty())
                out << head << " " << text << "\n";
            else
            {
                out << head << "\n";
                item.Dump(out, indent + 2);
            }
        }
    }

    std::ostream& operator<<(std::ostream& out, const ASDFnode& node)
    {
        node.Dump(out);
        return out;
    }
#pragma endregion
}
