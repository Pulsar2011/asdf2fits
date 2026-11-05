//
//  ASDFhandle.h
//  asdf2fits
//
//  Private: RAII wrappers around libasdf handles and the column layout shared by the table code.
//
//  Written by Claude (Anthropic's AI model) for William Gillard.
//  Copyright (c) 2026 William Gillard.
//
//  SPDX-License-Identifier: BSD-3-Clause
//  Licence text: see the LICENSE file at the root of this project.

/**
 * @file ASDFhandle.h
 * @brief RAII wrappers around libasdf handles and the column layout shared by the table code.
 * @note Private header: it isn't installed with the library.
 */
#ifndef _DSL_ASDFhandle_
#define _DSL_ASDFhandle_

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include <asdf.h>

#include <DSTfits/FITSexception.h>

/**
 * @namespace DSL::asdfdetail
 * @brief Implementation details of asdf2fits: conversion between ASDF and DeepSkyFits objects.
 */
namespace DSL::asdfdetail
{
    /** @brief Deleter of FilePtr. */
    struct FileCloser    { /** @brief Close a file. @param f File to close. */ void operator()(asdf_file_t* f)    const { asdf_close(f); } };
    /** @brief Deleter of ValuePtr. */
    struct ValueDeleter  { /** @brief Destroy a value. @param v Value to destroy. */ void operator()(asdf_value_t* v)   const { asdf_value_destroy(v); } };
    /** @brief Deleter of NdarrayPtr. */
    struct NdarrayDeleter{ /** @brief Destroy an ndarray. @param n Ndarray to destroy. */ void operator()(asdf_ndarray_t* n) const { asdf_ndarray_destroy(n); } };

    /** @brief Owning pointer to a libasdf file, closed on destruction. */
    using FilePtr    = std::unique_ptr<asdf_file_t, FileCloser>;
    /** @brief Owning pointer to a libasdf value, destroyed on destruction. */
    using ValuePtr   = std::unique_ptr<asdf_value_t, ValueDeleter>;
    /** @brief Owning pointer to a libasdf ndarray, destroyed on destruction. */
    using NdarrayPtr = std::unique_ptr<asdf_ndarray_t, NdarrayDeleter>;

    /**
     * @brief Text of the last libasdf error on a file.
     * @param file File to query; null when no file could be opened.
     * @return The error text, or a generic message if libasdf reports none.
     */
    inline std::string lastError(asdf_file_t* file)
    {
        const char* text = asdf_error(file);
        return text ? std::string(text) : std::string("unknown libasdf error");
    }

    /**
     * @brief Tell whether a value carries a tag containing the given text.
     * @param value Value of the ASDF tree.
     * @param text  Text to look for anywhere in the tag.
     * @return True if the value has a tag and the tag contains @p text.
     */
    inline bool hasTag(asdf_value_t* value, const char* text)
    {
        const char* tag = asdf_value_tag(value);
        return tag && std::string(tag).find(text) != std::string::npos;
    }

    /** @brief Layout of one table column, inside a row or in its own ndarray. */
    struct ColumnSpec
    {
        std::string name;               //!< Column name
        std::string unit;               //!< Physical unit; empty if none
        asdf_scalar_datatype_t type = ASDF_DATATYPE_UNKNOWN;    //!< Stored element type
        bool     swap      = false;     //!< Stored in the other byte order than this machine
        uint64_t elemBytes = 0;         //!< One element (for a string: all its characters)
        uint64_t repeat    = 1;         //!< Elements per row
        uint64_t offset    = 0;         //!< Start of the column inside a structured row
        asdf_ndarray_t* source = nullptr;   //!< ndarray holding the bytes; owned by the caller
    };
}

#endif
