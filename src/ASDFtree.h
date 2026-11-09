//
//  ASDFtree.h
//  asdf2fits
//
//  Private: conversion between libasdf values and ASDFnode trees, and writing of a tree as YAML.
//
//  Written by Claude (Anthropic's AI model) for William Gillard.
//  Copyright (c) 2026 William Gillard.
//
//  SPDX-License-Identifier: BSD-3-Clause
//  Licence text: see the LICENSE file at the root of this project.

/**
 * @file ASDFtree.h
 * @brief Conversion of libasdf values to ASDFnode trees, and writing of a tree as an ASDF YAML document.
 * @note Private header: it isn't installed with the library.
 */
#ifndef _DSL_ASDFtree_
#define _DSL_ASDFtree_

#include <functional>
#include <memory>
#include <string>

#include "ASDFnode.h"
#include "ASDFsource.h"

namespace DSL::asdfdetail
{
    /** @brief Tag of the root of an ASDF tree. */
    inline constexpr const char* kRootTag = "tag:stsci.edu:asdf/core/asdf-1.1.0";
    /** @brief Tag of an ndarray. */
    inline constexpr const char* kNdarrayTag = "tag:stsci.edu:asdf/core/ndarray-1.1.0";

    /**
     * @brief Convert a value of an open file, and everything below it, to a node.
     * @details Tags are kept, except the YAML core ones (tag:yaml.org,2002:...), which the kind of the node
     * already tells. Aliases are resolved: an anchored value appears at each of its aliases. Ndarrays become
     * handles on @p file, which describe them without reading their data.
     * @param value The value; its file is locked by the caller.
     * @param file  The open file, kept by the ndarray handles; may be null when the handles won't be read.
     * @param path  Path of @p value in the tree.
     * @return The node.
     */
    ASDFnode readTree(asdf_value_t* value, const std::shared_ptr<FileSource>& file, const std::string& path);

    /**
     * @brief Called for each ndarray node written by writeDocument().
     * @details Returns the node written in its place: an ndarray description, tagged as such.
     */
    using ArrayWriter = std::function<ASDFnode(const std::string& path, const ASDFnode& node)>;

    /**
     * @brief Write a tree as an ASDF YAML document, in block style.
     * @details Tags under tag:stsci.edu:asdf/ are written with the "!" handle, others in full. Strings that
     * would read back as another type are quoted.
     * @param root   Root of the tree; written with the core/asdf tag if it has none.
     * @param arrays Gives the node written in place of each ndarray node.
     * @return The document, from the "#ASDF" line to the "..." end marker.
     */
    std::string writeDocument(const ASDFnode& root, const ArrayWriter& arrays);
}

#endif
