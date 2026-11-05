//
//  ASDFheader.h
//  asdf2fits
//
//  Private: conversion between the ASDF tree and FITShdu keywords.
//
//  Written by Claude (Anthropic's AI model) for William Gillard.
//  Copyright (c) 2026 William Gillard.
//
//  SPDX-License-Identifier: BSD-3-Clause
//  Licence text: see the LICENSE file at the root of this project.

#ifndef _DSL_ASDFheader_
#define _DSL_ASDFheader_

#include <string>
#include <vector>

#include <DSTfits/FITShdu.h>

#include "ASDFhandle.h"

namespace DSL::asdfdetail
{
    /** @brief True for a tree value that isn't turned into keywords because it holds or is data. */
    bool isDataValue(asdf_value_t* value);

    /** @brief True for GWCS objects and transforms: their arrays aren't HDUs, their scalars aren't keywords. */
    bool isOpaqueValue(asdf_value_t* value);

    /** @brief Read a list of [key, value, comment] cards. */
    void readCards(asdf_value_t* cards, FITShdu& hdu);

    /** @brief Store one scalar of the tree as a keyword; anything else is ignored. */
    void readScalar(asdf_value_t* value, const std::string& key, FITShdu& hdu);

    /** @brief Turn the scalars below a mapping into keywords named PREFIX.KEY; data are skipped. */
    void readTreeKeywords(asdf_value_t* mapping, const std::string& prefix, FITShdu& hdu);

    /** @brief Set BITPIX, BSCALE, BZERO and NAXISn as FITSimg sets them for this pixel type. */
    void setImageStructure(FITShdu& hdu, asdf_scalar_datatype_t type, const std::vector<uint64_t>& naxis);

    /** @brief Set the keywords describing a binary table. */
    void setTableStructure(FITShdu& hdu, const std::vector<ColumnSpec>& columns, uint64_t rows);

    /** @brief Remove BSCALE, BZERO and the NAXISn left over from another shape. */
    void clearImageStructure(FITShdu& hdu);

    /** @brief Build the card list of a header; the keys in \c first come first, in that order. */
    asdf_sequence_t* writeCards(asdf_file_t* file, const FITShdu& hdu, const std::vector<std::string>& first);
}

#endif
