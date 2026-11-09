//
//  ASDFsource.h
//  asdf2fits
//
//  Private: the open file shared by ASDFmanager and the ndarray handles, and what a handle knows of its array.
//
//  Written by Claude (Anthropic's AI model) for William Gillard.
//  Copyright (c) 2026 William Gillard.
//
//  SPDX-License-Identifier: BSD-3-Clause
//  Licence text: see the LICENSE file at the root of this project.

/**
 * @file ASDFsource.h
 * @brief The open file shared by ASDFmanager and the ASDFarray handles, and the description of an array in it.
 * @note Private header: it isn't installed with the library.
 */
#ifndef _DSL_ASDFsource_
#define _DSL_ASDFsource_

#include <memory>
#include <mutex>
#include <string>
#include <vector>

#include "ASDFarray.h"
#include "ASDFhandle.h"
#include "ASDFnode.h"

namespace DSL::asdfdetail
{
    /**
     * @brief An open ASDF file, shared by its ASDFmanager and every ndarray handle read from it.
     * @details The file stays open as long as one of them refers to it. libasdf isn't thread-safe: every use
     * of @c file takes @c mtx.
     */
    struct FileSource
    {
        std::mutex mtx;         //!< Taken around every libasdf call on the file
        FilePtr file;           //!< The open file
        std::string path;       //!< Path of the file, for messages
        std::string tmpDir;     //!< Directory for decompressing to disk; libasdf keeps a pointer to it
    };

    /** @brief Part of an ndarray, in the ASDF (numpy) order: the slowest varying axis first. */
    struct Tile
    {
        std::vector<uint64_t> origin;   //!< First element, all axes; empty for the whole array
        std::vector<uint64_t> size;     //!< Length along each axis; empty for the whole array
        size_t dropped = 0;             //!< Leading axes, of length 1, that the image leaves out

        /** @brief Tell whether the tile is the whole array. @return True if it is. */
        bool whole() const { return size.empty(); }
    };

    /** @brief What a handle knows of its array, as found in the file. Shared by the handles of the same array. */
    struct ArraySource
    {
        std::shared_ptr<FileSource> file;   //!< File holding the array
        std::string path;                   //!< Path of the ndarray value in the tree (an ASDFnode path)
        ASDFnode desc;                      //!< The ndarray value as written, with its tag; written back for a verbatim copy
        std::vector<uint64_t> shape;        //!< Shape, slowest axis first; empty if unknown
        asdf_scalar_datatype_t type = ASDF_DATATYPE_UNKNOWN; //!< Element type, if libasdf can read the array
        std::string datatype;               //!< Name of the element type, as returned by ASDFarray::DataType()
        uint64_t elementBytes = 0;          //!< Size of one element; 0 if unknown
        int64_t block = -1;                 //!< Index of its block in the file; -1 for inline or external data
    };

    /**
     * @brief Describe an ndarray value of an open file, without reading its data.
     * @param file  File holding the value.
     * @param value The ndarray value.
     * @param path  Its path in the tree.
     * @param desc  The value converted to a node.
     * @return The description.
     */
    std::shared_ptr<ArraySource> describeArray(const std::shared_ptr<FileSource>& file, asdf_value_t* value,
                                               const std::string& path, ASDFnode desc);

    /**
     * @brief Find the value at a path of a file's tree.
     * @param file An open file, locked by the caller.
     * @param path ASDFnode path of the value.
     * @return The value, or null if the path leads nowhere.
     */
    ValuePtr valueAt(asdf_file_t* file, const std::string& path);

    /** @brief Access to the private state of ASDFarray, for the library. */
    struct ArrayAccess
    {
        /**
         * @brief Make a handle on an array of a file.
         * @param source Description of the array.
         * @return The handle, with no view and nothing loaded.
         */
        static std::shared_ptr<ASDFarray> make(std::shared_ptr<const ArraySource> source);

        /** @brief Get where a handle's array is. @param array The handle. @return The description, or null for a handle in memory. */
        static const ArraySource* source(const ASDFarray& array) { return array.fsource.get(); }

        /** @brief Get the shared description of a handle's array. @param array The handle. @return The description, or null. */
        static std::shared_ptr<const ArraySource> sharedSource(const ASDFarray& array) { return array.fsource; }

        /** @brief Get the pixels in memory. @param array The handle. @return The image, or null if not loaded. */
        static const std::shared_ptr<FITScube>& cube(const ASDFarray& array) { return array.fcube; }

        /**
         * @brief Get the view of a handle in its file.
         * @param array The handle.
         * @return The tile it shows, all axes, in the ASDF order.
         */
        static Tile tile(const ASDFarray& array);

        /**
         * @brief Get the ndarray of a handle's array.
         * @param source Where the array is; its file must be locked by the caller.
         * @return The ndarray.
         * @throws FITSexception READ_ERROR if the array can't be found or read.
         */
        static NdarrayPtr open(const ArraySource& source);
    };
}

#endif
