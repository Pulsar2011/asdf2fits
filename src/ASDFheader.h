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

/**
 * @file ASDFheader.h
 * @brief Conversion between the ASDF tree and FITShdu keywords.
 * @note Private header: it isn't installed with the library.
 */
#ifndef _DSL_ASDFheader_
#define _DSL_ASDFheader_

#include <string>
#include <vector>

#include <DSTfits/FITShdu.h>

#include "ASDFhandle.h"

namespace DSL::asdfdetail
{
    /**
     * @brief Tell whether a tree value is data rather than keywords.
     * @param value Value of the ASDF tree.
     * @return True for a value tagged as an ndarray, a table or an astropy FITS object.
     */
    bool isDataValue(asdf_value_t* value);

    /**
     * @brief Tell whether a tree value is a GWCS object or a transform.
     * @details These values are skipped: their arrays aren't HDUs and their scalars aren't keywords.
     * @param value Value of the ASDF tree.
     * @return True if the tag of the value contains "gwcs/" or "/transform".
     */
    bool isOpaqueValue(asdf_value_t* value);

    /**
     * @brief Read a list of [key, value, comment] cards into a header.
     * @details The value and the comment of a card are optional. COMMENT and HISTORY cards are
     * joined into one entry, as FITShdu does for a FITS file, and text under a blank keyword counts
     * as COMMENT. Cards without a text key, such as the [] astropy writes for a blank card, are skipped.
     * @param cards Sequence of cards; anything else is ignored.
     * @param hdu   Header receiving the keywords.
     */
    void readCards(asdf_value_t* cards, FITShdu& hdu);

    /**
     * @brief Store one scalar of the tree as a keyword.
     * @param value Value to store; a mapping, a sequence or a null value is ignored.
     * @param key   Name of the keyword, used as it is.
     * @param hdu   Header receiving the keyword.
     */
    void readScalar(asdf_value_t* value, const std::string& key, FITShdu& hdu);

    /**
     * @brief Turn the scalars below a mapping into keywords named PREFIX.KEY.
     * @details Nested mappings extend the prefix, and keys are put in upper case. Data, GWCS objects,
     * transforms and sequences are skipped.
     * @param mapping Mapping to read; anything else is ignored.
     * @param prefix  Start of the keyword names, as it is; empty for none.
     * @param hdu     Header receiving the keywords.
     */
    void readTreeKeywords(asdf_value_t* mapping, const std::string& prefix, FITShdu& hdu);

    /**
     * @brief Set BITPIX, BSCALE, BZERO and NAXISn as FITSimg sets them for this pixel type.
     * @details Any previous BSCALE, BZERO and NAXISn are removed first. A type with no image form
     * gets BITPIX = 8, without BSCALE or BZERO.
     * @param hdu   Header to update.
     * @param type  Pixel type.
     * @param naxis Axis lengths in FITS order: NAXIS1 first.
     */
    void setImageStructure(FITShdu& hdu, asdf_scalar_datatype_t type, const std::vector<uint64_t>& naxis);

    /**
     * @brief Set the keywords describing a binary table.
     * @details Sets XTENSION, BITPIX, NAXIS, NAXIS1, NAXIS2, PCOUNT, GCOUNT, TFIELDS and, per column,
     * TTYPEn, TFORMn and TUNITn if it has a unit. Removes any BSCALE, BZERO, BLANK and other NAXISn.
     * @param hdu     Header to update.
     * @param columns Layout of the columns, in table order.
     * @param rows    Number of rows.
     */
    void setTableStructure(FITShdu& hdu, const std::vector<ColumnSpec>& columns, uint64_t rows);

    /**
     * @brief Remove BSCALE, BZERO, NAXIS and the NAXISn left over from another shape.
     * @param hdu Header to update.
     */
    void clearImageStructure(FITShdu& hdu);

    /**
     * @brief Build the card list of a header.
     * @details The keys in @p first come first, in that order, followed by the other keywords.
     * COMMENT and HISTORY come last, split into cards of 72 characters.
     * @param file  File the list is created in.
     * @param hdu   Header to write.
     * @param first Keys to put first; those missing from the header are skipped.
     * @return New sequence of [key, value, comment] cards, without a comment if the keyword has none.
     */
    asdf_sequence_t* writeCards(asdf_file_t* file, const FITShdu& hdu, const std::vector<std::string>& first);
}

#endif
