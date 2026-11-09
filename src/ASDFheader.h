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

#include <map>
#include <string>
#include <vector>

#include <DSTfits/FITShdu.h>

#include "ASDFhandle.h"
#include "ASDFnode.h"

namespace DSL::asdfdetail
{
    /** @brief Path in the tree of the value behind each keyword of a header. */
    using KeyPaths = std::map<std::string, std::string>;

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

    /** @brief As isDataValue(), for a node. @param node Node of the tree. @return True for data. */
    bool isDataNode(const ASDFnode& node);

    /** @brief As isOpaqueValue(), for a node. @param node Node of the tree. @return True for a GWCS object or a transform. */
    bool isOpaqueNode(const ASDFnode& node);

    /**
     * @brief Tell whether a node holds, at any depth, something that becomes an HDU.
     * @param node Node of the tree.
     * @return True if it is or contains data outside any GWCS object or transform.
     */
    bool holdsData(const ASDFnode& node);

    /**
     * @brief Read a list of [key, value, comment] cards into a header.
     * @details The value and the comment of a card are optional. COMMENT and HISTORY cards are
     * joined into one entry, as FITShdu does for a FITS file, and text under a blank keyword counts
     * as COMMENT. Cards without a text key, such as the [] astropy writes for a blank card, are skipped.
     * @param cards Sequence of cards; anything else is ignored.
     * @param hdu   Header receiving the keywords.
     */
    void readCards(const ASDFnode& cards, FITShdu& hdu);

    /**
     * @brief Store one scalar of the tree as a keyword.
     * @param value Value to store; a container, an ndarray or a null value is ignored.
     * @param key   Name of the keyword, used as it is.
     * @param hdu   Header receiving the keyword.
     * @return True if a keyword was set.
     */
    bool readScalar(const ASDFnode& value, const std::string& key, FITShdu& hdu);

    /**
     * @brief Turn the scalars below a mapping into keywords named PREFIX.KEY.
     * @details Nested mappings extend the prefix, and keys are put in upper case. Data, GWCS objects,
     * transforms and sequences are skipped.
     * @param mapping Mapping to read; anything else is ignored.
     * @param prefix  Start of the keyword names, as it is; empty for none.
     * @param path    Path of @p mapping in the tree, to fill @p paths.
     * @param hdu     Header receiving the keywords.
     * @param paths   Receives the path of the value of each keyword set; may be null.
     */
    void readTreeKeywords(const ASDFnode& mapping, const std::string& prefix, const std::string& path, FITShdu& hdu, KeyPaths* paths = nullptr);

    /**
     * @brief Turn the root of the tree into keywords.
     * @details Entries holding data, GWCS objects, transforms, asdf_library and history are skipped. A mapping
     * MAP gives keywords MAP.KEY, a scalar gives its own.
     * @param root  Root of the tree; ignored if it isn't a mapping.
     * @param hdu   Header receiving the keywords.
     * @param paths Receives the path of the value of each keyword set; may be null.
     */
    void readRootKeywords(const ASDFnode& root, FITShdu& hdu, KeyPaths* paths = nullptr);

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

    /**
     * @brief Tell whether a keyword describes the data rather than the content of an image.
     * @details SIMPLE, XTENSION, BITPIX, NAXIS, NAXISn, BSCALE, BZERO, BLANK, EXTEND, PCOUNT, GCOUNT, CHECKSUM
     * and DATASUM. Outside the astropy FITS layout, EXTNAME too, as it comes from the tree path.
     * @param key        Keyword.
     * @param withExtname Count EXTNAME as structural.
     * @return True for those keywords.
     */
    bool isStructuralKey(const std::string& key, bool withExtname);

    /** @brief How the header of an image maps to the tree it was read from. */
    struct HeaderLink
    {
        std::string cardsPath;      //!< Astropy FITS layout: path of the card list; empty for any other layout
        std::string keywordsPath;   //!< Other layouts: path of the mapping whose scalars are keywords
        bool hasKeywords = false;   //!< @c keywordsPath holds keywords (false for data at the root of the tree)
        bool rootKeywords = false;  //!< The scalars at the root of the tree are keywords too (the primary HDU)
        std::string sinkPath;       //!< Other layouts: mapping receiving the keywords that have no place in the tree
        KeyPaths keyPaths;          //!< Other layouts: path of the value of each keyword read from the tree
    };

    /**
     * @brief Copy into the tree the changes made to a header since a snapshot of it.
     * @details In the astropy FITS layout, the cards of changed keywords are replaced, those of removed
     * keywords are removed and new keywords are added at the end; COMMENT and HISTORY are written again if
     * they changed. In other layouts, a changed keyword updates the value it was read from, keeping its tag,
     * a removed keyword removes it, and a new keyword is added to the sink mapping (as nested mappings for
     * a dotted name, in lower case). Structural keywords are left out (see isStructuralKey()).
     * @param tree     Tree to update.
     * @param link     How the header maps to @p tree; @c keyPaths receives the new keywords.
     * @param hdu      The header now.
     * @param snapshot The header when it last matched @p tree.
     */
    void syncHeader(ASDFnode& tree, HeaderLink& link, const FITShdu& hdu, const FITShdu& snapshot);

    /**
     * @brief Remove the cards of a keyword from an astropy card list.
     * @param cards Card list.
     * @param key   Keyword; for COMMENT, blank-key cards holding text go too.
     */
    void removeCards(ASDFnode& cards, const std::string& key);

    /**
     * @brief Read again from the tree the keywords of a header.
     * @details Keywords that are no longer in the tree are removed; structural keywords are left as they are.
     * @param tree Tree to read.
     * @param link How the header maps to @p tree; @c keyPaths is filled again.
     * @param hdu  Header to update.
     */
    void rebuildHeader(const ASDFnode& tree, HeaderLink& link, FITShdu& hdu);
}

#endif
