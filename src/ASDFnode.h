//
//  ASDFnode.h
//  asdf2fits
//
//  Written by Claude (Anthropic's AI model) for William Gillard.
//  Copyright (c) 2026 William Gillard.
//
//  SPDX-License-Identifier: BSD-3-Clause
//  Licence text: see the LICENSE file at the root of this project.

/**
 * @file ASDFnode.h
 * @brief A value of an ASDF tree: scalars, sequences and mappings with their YAML tags, and ndarrays.
 */
#ifndef _DSL_ASDFnode_
#define _DSL_ASDFnode_

#include <concepts>
#include <cstdint>
#include <functional>
#include <iosfwd>
#include <limits>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace DSL
{
    class ASDFarray;

    /**
     * @class DSL::ASDFnode ASDFnode.h "asdf2fits/ASDFnode.h"
     * @brief One value of an ASDF tree, with everything below it.
     * @details A node is null, a scalar (boolean, integer, real or string), a sequence, a mapping or an ndarray.
     * Every node keeps its YAML tag, so values that asdf2fits doesn't interpret (GWCS objects, units, times,
     * transforms) are written back unchanged. Mappings keep their keys in file order.
     *
     * An ndarray node holds a handle to an ASDFarray, which describes the array without reading it. Copying a
     * node copies the tree and shares these handles; Clone() also copies the handles, so that the copy is
     * fully independent.
     *
     * Paths name a node below another one: mapping keys and sequence indices joined by '/', as in
     * "roman/meta/exposure/start_time" or "fits/0/header". A '/' inside a key is written "~1" and a '~' is
     * written "~0", as in a JSON pointer. The empty path names the node itself.
     *
     * Failures throw FITSexception.
     */
    class ASDFnode
    {
    public:
        /** @brief What a node holds. */
        enum class Kind
        {
            null,               //!< No value
            boolean,            //!< true or false
            integer,            //!< Signed integer (int64_t), or an unsigned one that libasdf reads as a small type
            unsigned_integer,   //!< Unsigned integer (uint64_t)
            real,               //!< Floating-point number
            string,             //!< Text, including tagged scalars such as times and units
            sequence,           //!< Ordered list of nodes
            mapping,            //!< Nodes named by keys, in file order
            array               //!< An ndarray, through an ASDFarray handle
        };

#pragma region * ctor
        /** @brief Make a null node. */
        ASDFnode() = default;

        /** @brief Make a null node. */
        ASDFnode(std::nullptr_t) {}

        /**
         * @brief Make a boolean.
         * @param value The value.
         */
        ASDFnode(bool value): fkind(Kind::boolean), fbool(value) {}

        /**
         * @brief Make an integer.
         * @param value The value.
         */
        template<std::signed_integral I>
        ASDFnode(I value): fkind(Kind::integer), fint(static_cast<int64_t>(value)) {}

        /**
         * @brief Make an unsigned integer.
         * @param value The value.
         */
        template<std::unsigned_integral U> requires (!std::same_as<U, bool>)
        ASDFnode(U value): fkind(Kind::unsigned_integer), fuint(static_cast<uint64_t>(value)) {}

        /**
         * @brief Make a real.
         * @param value The value.
         */
        template<std::floating_point F>
        ASDFnode(F value): fkind(Kind::real), freal(static_cast<double>(value)) {}

        /**
         * @brief Make a string.
         * @param text The text; a null pointer makes a null node.
         */
        ASDFnode(const char* text);

        /**
         * @brief Make a string.
         * @param text The text.
         */
        ASDFnode(std::string text): fkind(Kind::string), ftext(std::move(text)) {}

        /**
         * @brief Make an ndarray node.
         * @param array Handle to the array; a null handle makes a null node.
         */
        ASDFnode(std::shared_ptr<ASDFarray> array);

        /**
         * @brief Make an empty mapping.
         * @param tag YAML tag of the mapping, in full; empty for none.
         * @return The mapping.
         */
        static ASDFnode Mapping(std::string tag = {});

        /**
         * @brief Make an empty sequence.
         * @param tag YAML tag of the sequence, in full; empty for none.
         * @return The sequence.
         */
        static ASDFnode Sequence(std::string tag = {});

#pragma endregion
#pragma region * Kind and tag
        /** @brief What the node holds. @return The kind. */
        Kind GetKind() const { return fkind; }

        bool IsNull()     const { return fkind == Kind::null; }                                         //!< @brief Tell whether the node is null. @return True if null.
        bool IsBool()     const { return fkind == Kind::boolean; }                                      //!< @brief Tell whether the node is a boolean. @return True if boolean.
        bool IsInteger()  const { return fkind == Kind::integer || fkind == Kind::unsigned_integer; }   //!< @brief Tell whether the node is an integer. @return True if integer.
        bool IsReal()     const { return fkind == Kind::real; }                                         //!< @brief Tell whether the node is a real. @return True if real.
        bool IsNumber()   const { return IsInteger() || IsReal(); }                                     //!< @brief Tell whether the node is a number. @return True if integer or real.
        bool IsString()   const { return fkind == Kind::string; }                                       //!< @brief Tell whether the node is a string. @return True if string.
        bool IsScalar()   const { return fkind != Kind::sequence && fkind != Kind::mapping && fkind != Kind::array; }  //!< @brief Tell whether the node is null or a scalar. @return True if neither a container nor an ndarray.
        bool IsSequence() const { return fkind == Kind::sequence; }                                     //!< @brief Tell whether the node is a sequence. @return True if sequence.
        bool IsMapping()  const { return fkind == Kind::mapping; }                                      //!< @brief Tell whether the node is a mapping. @return True if mapping.
        bool IsArray()    const { return fkind == Kind::array; }                                        //!< @brief Tell whether the node is an ndarray. @return True if ndarray.

        /** @brief Get the YAML tag. @return The tag in full, such as "tag:stsci.edu:asdf/time/time-1.4.0"; empty if none. */
        const std::string& Tag() const { return ftag; }

        /** @brief Set the YAML tag. @param tag The tag in full; empty to remove it. */
        void SetTag(std::string tag) { ftag = std::move(tag); }

        /**
         * @brief Tell whether the tag contains a text.
         * @param text Text to look for, such as "gwcs/".
         * @return True if the node has a tag containing @p text.
         */
        bool HasTag(const std::string& text) const { return !ftag.empty() && ftag.find(text) != std::string::npos; }

#pragma endregion
#pragma region * Scalars
        /**
         * @brief Get the value of a scalar.
         * @details T can be bool, any integer type, float, double or std::string. Integers convert to reals.
         * @tparam T Type to return.
         * @return The value.
         * @throws FITSexception BAD_DATATYPE if the node holds another kind of value, NUM_OVERFLOW if an
         * integer doesn't fit in @p T.
         */
        template<typename T>
        T Get() const;

        /**
         * @brief Get a scalar as text.
         * @return The value as it is written in YAML (".nan", "true", "42"...); empty for a null node or a container.
         */
        std::string Text() const;

#pragma endregion
#pragma region * ndarray
        /** @brief Get the ndarray handle. @return The handle, or null if the node isn't an ndarray. */
        const std::shared_ptr<ASDFarray>& Array() const { return farray; }

#pragma endregion
#pragma region * Containers
        /** @brief Get the number of items. @return The items of a sequence or the entries of a mapping; 0 otherwise. */
        size_t Size() const { return fitems.size(); }

        /**
         * @brief Tell whether a mapping has a key.
         * @param key Key to look for.
         * @return True if the node is a mapping with this key.
         */
        bool Contains(const std::string& key) const { return Find(key) != nullptr; }

        /**
         * @brief Find the value of a key.
         * @param key Key to look for.
         * @return The value, or null if the node isn't a mapping or has no such key.
         */
        const ASDFnode* Find(const std::string& key) const;
        ASDFnode*       Find(const std::string& key);       //!< @copydoc Find(const std::string&) const

        /**
         * @brief Get the value of a key, adding a null value if it is missing.
         * @details A null node becomes an empty mapping first.
         * @param key Key of the value.
         * @return The value.
         * @throws FITSexception BAD_DATATYPE if the node is neither a mapping nor null.
         */
        ASDFnode& operator[](const std::string& key);

        /**
         * @brief Get the value of a key.
         * @param key Key of the value.
         * @return The value.
         * @throws FITSexception KEY_NO_EXIST if the node isn't a mapping or has no such key.
         */
        const ASDFnode& operator[](const std::string& key) const;

        /**
         * @brief Get an item of a sequence, or the value of the n-th entry of a mapping.
         * @param index Index of the item, from 0.
         * @return The item.
         * @throws FITSexception BAD_DIMEN if there is no such item.
         */
        ASDFnode&       operator[](size_t index);
        const ASDFnode& operator[](size_t index) const;     //!< @copydoc operator[](size_t)

        /**
         * @brief Get the keys of a mapping.
         * @return The keys in file order, matching Items(); empty for any other node.
         */
        const std::vector<std::string>& Keys() const { return fkeys; }

        /**
         * @brief Get the items of a sequence, or the values of a mapping.
         * @return The items, in order.
         */
        const std::vector<ASDFnode>& Items() const { return fitems; }
        std::vector<ASDFnode>&       Items()       { return fitems; }   //!< @copydoc Items() const

        /**
         * @brief Add an item at the end of a sequence.
         * @details A null node becomes an empty sequence first.
         * @param value Item to add.
         * @return The added item.
         * @throws FITSexception BAD_DATATYPE if the node is neither a sequence nor null.
         */
        ASDFnode& Append(ASDFnode value);

        /**
         * @brief Remove a key from a mapping.
         * @param key Key to remove.
         * @return True if it was there.
         */
        bool Erase(const std::string& key);

        /**
         * @brief Remove an item from a sequence or a mapping.
         * @param index Index of the item, from 0.
         * @return True if it was there.
         */
        bool EraseAt(size_t index);

#pragma endregion
#pragma region * Paths
        /**
         * @brief Find the node at a path.
         * @param path Path below this node; empty for the node itself.
         * @return The node, or null if the path leads nowhere.
         */
        const ASDFnode* FindPath(const std::string& path) const;
        ASDFnode*       FindPath(const std::string& path);      //!< @copydoc FindPath(const std::string&) const

        /**
         * @brief Get the node at a path.
         * @param path Path below this node.
         * @return The node.
         * @throws FITSexception KEY_NO_EXIST if the path leads nowhere.
         */
        const ASDFnode& At(const std::string& path) const;
        ASDFnode&       At(const std::string& path);            //!< @copydoc At(const std::string&) const

        /**
         * @brief Set the node at a path, adding the missing mappings on the way.
         * @details A sequence index on the way must exist. The tag of a scalar replaced by a scalar is kept, so
         * that a time or a unit stays one.
         * @param path  Path below this node; not empty.
         * @param value Value to set.
         * @return The node set.
         * @throws FITSexception BAD_DATATYPE if a node on the way is a scalar, KEY_NO_EXIST if a sequence index is out of range.
         */
        ASDFnode& SetPath(const std::string& path, ASDFnode value);

        /**
         * @brief Remove the node at a path from its parent.
         * @param path Path below this node; not empty.
         * @return True if it was there.
         */
        bool ErasePath(const std::string& path);

        /**
         * @brief Join path parts, escaping '~' and '/'.
         * @param parts Keys and indices.
         * @return The path.
         */
        static std::string JoinPath(const std::vector<std::string>& parts);

        /**
         * @brief Add a part at the end of a path.
         * @param path Path to extend; may be empty.
         * @param part Key or index to add, unescaped.
         * @return The longer path.
         */
        static std::string ChildPath(const std::string& path, const std::string& part);

        /**
         * @brief Split a path into its parts, unescaping '~0' and '~1'.
         * @param path Path to split.
         * @return Its keys and indices; empty for the empty path.
         */
        static std::vector<std::string> SplitPath(const std::string& path);

        /**
         * @brief Get the path of the parent.
         * @param path Path of a node; not empty.
         * @return The path without its last part.
         */
        static std::string ParentPath(const std::string& path);

#pragma endregion
#pragma region * Tree
        /**
         * @brief Visit this node and every node below it, parents first.
         * @param visitor Called with the path of each node below this one (empty for this node) and the node;
         * returns false to skip what is below that node. Don't add or remove items of the node's parent.
         */
        void Visit(const std::function<bool(const std::string& path, ASDFnode& node)>& visitor);
        void Visit(const std::function<bool(const std::string& path, const ASDFnode& node)>& visitor) const;   //!< @copydoc Visit

        /**
         * @brief Make a fully independent copy.
         * @details Unlike the copy constructor, ndarray handles are copied too (see ASDFarray::Clone()).
         * @return The copy.
         */
        ASDFnode Clone() const;

        /**
         * @brief Compare two trees: kinds, values, tags, keys and items; ndarrays by their handle.
         * @details Integers are equal if their values are, whether signed or unsigned.
         * @param other Tree to compare with.
         * @return True if they are the same.
         */
        bool operator==(const ASDFnode& other) const;

        /**
         * @brief Print the tree in a YAML-like form.
         * @param out    Stream to print to.
         * @param indent Indentation of the first level.
         */
        void Dump(std::ostream& out, int indent = 0) const;

#pragma endregion
    private:
        [[noreturn]] void typeError(const char* wanted) const;  //!< Throw the exception of Get() for a mismatch
        [[noreturn]] void rangeError(const char* wanted) const; //!< Throw the exception of Get() for an overflow

        Kind fkind = Kind::null;            //!< What the node holds
        std::string ftag;                   //!< YAML tag in full; empty if none
        bool     fbool = false;             //!< Value of a boolean
        int64_t  fint  = 0;                 //!< Value of an integer
        uint64_t fuint = 0;                 //!< Value of an unsigned_integer
        double   freal = 0.;                //!< Value of a real
        std::string ftext;                  //!< Value of a string
        std::vector<std::string> fkeys;     //!< Keys of a mapping, matching fitems
        std::vector<ASDFnode> fitems;       //!< Items of a sequence or values of a mapping
        std::shared_ptr<ASDFarray> farray;  //!< Handle of an ndarray
    };

    template<typename T>
    T ASDFnode::Get() const
    {
        if constexpr (std::same_as<T, bool>)
        {
            if(fkind != Kind::boolean)
                typeError("a boolean");
            return fbool;
        }
        else if constexpr (std::integral<T>)
        {
            if(fkind == Kind::integer)
            {
                if(!std::in_range<T>(fint))
                    rangeError("this integer type");
                return static_cast<T>(fint);
            }
            if(fkind == Kind::unsigned_integer)
            {
                if(!std::in_range<T>(fuint))
                    rangeError("this integer type");
                return static_cast<T>(fuint);
            }
            typeError("an integer");
        }
        else if constexpr (std::floating_point<T>)
        {
            switch(fkind)
            {
                case Kind::real:             return static_cast<T>(freal);
                case Kind::integer:          return static_cast<T>(fint);
                case Kind::unsigned_integer: return static_cast<T>(fuint);
                default:                     typeError("a number");
            }
        }
        else if constexpr (std::same_as<T, std::string>)
        {
            if(fkind != Kind::string)
                typeError("a string");
            return ftext;
        }
        else
            static_assert(!sizeof(T), "ASDFnode::Get supports bool, integers, reals and std::string");
    }

    /**
     * @brief Print a tree in a YAML-like form.
     * @param out  Stream to print to.
     * @param node Tree to print.
     * @return The stream.
     */
    std::ostream& operator<<(std::ostream& out, const ASDFnode& node);
}

#endif
