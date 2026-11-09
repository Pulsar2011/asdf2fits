//
//  ASDFwrite.h
//  asdf2fits
//
//  Private: writing of an ASDFnode tree, with its arrays, as an ASDF file.
//
//  Written by Claude (Anthropic's AI model) for William Gillard.
//  Copyright (c) 2026 William Gillard.
//
//  SPDX-License-Identifier: BSD-3-Clause
//  Licence text: see the LICENSE file at the root of this project.

/**
 * @file ASDFwrite.h
 * @brief Writing of an ASDFnode tree, with its arrays, as an ASDF file.
 * @note Private header: it isn't installed with the library.
 */
#ifndef _DSL_ASDFwrite_
#define _DSL_ASDFwrite_

#include <string>

#include "ASDFimage.h"
#include "ASDFnode.h"

namespace DSL::asdfdetail
{
    /**
     * @brief Check a compression option.
     * @param compression "keep", "", "zlib", "bzp2" or "lz4".
     * @throws FITSexception BAD_OPTION for anything else.
     */
    void checkCompression(const std::string& compression);

    /**
     * @brief Make a history entry dated now, signed by asdf2fits.
     * @param description What was done.
     * @return The entry, tagged as a core/history_entry.
     */
    ASDFnode historyEntry(const std::string& description);

    /**
     * @brief Add an entry to the history of a tree.
     * @details Goes to history/entries, made if needed; to history itself if it is an old-style sequence.
     * @param tree  Root of the tree.
     * @param entry Entry to add.
     */
    void appendHistory(ASDFnode& tree, ASDFnode entry);

    /**
     * @brief Write a tree as an ASDF file.
     * @details Each ndarray node becomes a block of the file:
     *  - the image at @p imagePath, from @p image;
     *  - an array loaded or set, from memory;
     *  - an array with a view, by reading its view;
     *  - any other array of a file, by copying its block as it is, without decompressing it; arrays sharing a
     *    block keep sharing it.
     * Inline and external arrays are written as they were. The file is written under a temporary name in the
     * same directory, then renamed.
     * @param tree        Tree to write.
     * @param imagePath   Path of the image's ndarray; empty for none.
     * @param image       Pixels of the image; those it borrows must stay valid until the end of the call.
     * @param compression Compression of the arrays written from memory: "keep" for their original one, "" for
     * none, or a compressor.
     * @param path        File to create.
     * @param replace     Overwrite an existing file.
     * @throws FITSexception FILE_NOT_CREATED if the file exists and @p replace is false, or if it can't be
     * written; BAD_OPTION for an unknown compression; READ_ERROR if an array can't be read; WRITE_ERROR if
     * libasdf refuses the tree or a block.
     */
    void writeTreeFile(const ASDFnode& tree, const std::string& imagePath, PixelBlock image,
                       const std::string& compression, const std::string& path, bool replace);

    /**
     * @brief Write a libasdf file, replacing any file at the path.
     * @details libasdf writes over an existing file without truncating it, which would leave the end of a
     * longer file, and its block index, behind. The file is written under a new temporary name in the same
     * directory, then renamed, which also leaves the old file whole if writing fails.
     * @param file File to write; the blocks it borrows must stay valid until the end of the call.
     * @param path File to create or replace.
     * @throws FITSexception FILE_NOT_CREATED if the file can't be written or renamed.
     */
    void writeReplacing(asdf_file_t* file, const std::string& path);
}

#endif
