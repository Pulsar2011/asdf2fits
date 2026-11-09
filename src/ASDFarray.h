//
//  ASDFarray.h
//  asdf2fits
//
//  Written by Claude (Anthropic's AI model) for William Gillard.
//  Copyright (c) 2026 William Gillard.
//
//  SPDX-License-Identifier: BSD-3-Clause
//  Licence text: see the LICENSE file at the root of this project.

/**
 * @file ASDFarray.h
 * @brief Lazy handle to an ndarray of an ASDF file.
 */
#ifndef _DSL_ASDFarray_
#define _DSL_ASDFarray_

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace DSL
{
    class FITScube;

    namespace asdfdetail
    {
        struct ArraySource;
        struct ArrayAccess;
    }

    /**
     * @class DSL::ASDFarray ASDFarray.h "asdf2fits/ASDFarray.h"
     * @brief Handle to an ndarray, read from its file only when asked.
     * @details The shape, the data type and the size of the array are known without reading it. Load() reads
     * the pixels as a FITSimg\<T\> and keeps them; Read() reads one tile and keeps nothing; Release() frees the
     * pixels kept by Load(). Nothing else reads the pixels, so arrays that are never asked for cost no memory.
     *
     * A handle can show only part of its array, a view made by Window() or Layer(): the view costs nothing
     * until it is loaded, and then only its pixels are read. Shapes, positions and sizes are given in the
     * FITS order, NAXIS1 (the fastest varying axis) first, and relative to the view.
     *
     * A handle read from a file keeps that file open, even after its ASDFmanager is closed. A handle made by
     * FromCube() lives in memory only.
     *
     * Pixel types are those of ASDFmanager::GetImage(): booleans load as uint8_t and 16-bit reals as float.
     * Pixels that are NaN are masked.
     *
     * A handle isn't meant to be used from several threads at once; reads from the same file are serialised.
     * Failures throw FITSexception.
     */
    class ASDFarray
    {
    public:
        ~ASDFarray();                                           //!< Release the pixels and the file
        ASDFarray(const ASDFarray&) = delete;                   //!< Use Clone()
        ASDFarray& operator=(const ASDFarray&) = delete;        //!< Use Clone()

#pragma region * Description, without reading
        /**
         * @brief Get the axis lengths.
         * @return The lengths of the view, NAXIS1 first.
         */
        std::vector<size_t> Axes() const;

        /**
         * @brief Get the shape, in the ASDF (numpy) order.
         * @return The lengths of the view, the slowest varying axis first.
         */
        std::vector<uint64_t> Shape() const;

        /** @brief Get the number of elements. @return The number of elements of the view. */
        uint64_t Elements() const;

        /** @brief Get the size of the data. @return The size in bytes of the view, uncompressed, in its file type. */
        uint64_t Bytes() const;

        /**
         * @brief Get the type of the elements, as written in the file.
         * @return "float32", "uint16", "bool8"..., "ascii" or "ucs4" for strings, "structured" for records.
         */
        std::string DataType() const;

        /** @brief Tell whether the array loads as an image. @return True for integers, reals and booleans. */
        bool IsImage() const;

        /**
         * @brief Get the compression of the array in its file.
         * @details Reads the header of its block the first time, not the data.
         * @return "zlib", "bzp2", "lz4", or empty if the data aren't compressed or the handle is in memory.
         */
        std::string Compression() const;

        /** @brief Tell whether the array comes from a file. @return False for a handle made by FromCube(). */
        bool IsInFile() const;

        /** @brief Tell whether the pixels are in memory. @return True after Load(), until Release(). */
        bool IsLoaded() const;

        /** @brief Tell whether the handle shows part of its array. @return True after Window() or Layer(). */
        bool HasView() const;

#pragma endregion
#pragma region * Data
        /**
         * @brief Read the pixels of the view, once, and keep them.
         * @details Changes made to the returned image are those written by ASDFcube::WriteASDF().
         * @return The image; the same one until Release().
         * @throws FITSexception BAD_DATATYPE if the array has no image form, READ_ERROR if reading fails.
         */
        std::shared_ptr<FITScube> Load();

        /**
         * @brief Read a tile of the view, and keep nothing.
         * @details Uses the pixels in memory if they are loaded.
         * @param origin First pixel of the tile, NAXIS1 first; missing axes start at 0.
         * @param size   Length of the tile along each axis, NAXIS1 first; missing or 0 lengths go to the end of the axis.
         * @return The tile, as an image with as many axes as the view.
         * @throws FITSexception BAD_DIMEN if the tile isn't inside the view, or as Load().
         */
        std::shared_ptr<FITScube> Read(const std::vector<size_t>& origin, const std::vector<size_t>& size = {}) const;

        /**
         * @brief Read one element.
         * @param pixel Position of the element, NAXIS1 first.
         * @return The value, converted to double.
         * @throws FITSexception BAD_DIMEN if the position isn't inside the view, or as Load().
         */
        double ValueAt(const std::vector<size_t>& pixel) const;

        /**
         * @brief Free the pixels kept by Load().
         * @details Changes made to them are lost: the next Load() reads the file again. Pixels that can't be
         * read again, those of a handle made by FromCube() or given to Set(), are kept.
         * @return True if the pixels were freed or weren't loaded; false if they were kept.
         */
        bool Release();

        /**
         * @brief Replace the content of the array.
         * @details The handle then lives in memory: it no longer reads its file.
         * @param cube New content; not null.
         * @throws FITSexception SHARED_NULPTR if @p cube is null.
         */
        void Set(std::shared_ptr<FITScube> cube);

#pragma endregion
#pragma region * Views and copies
        /**
         * @brief Make a handle showing part of this one.
         * @details Nothing is read: the new handle reads only its window when it is loaded. If this handle is
         * loaded, the window is copied from memory instead.
         * @param origin First pixel of the window, NAXIS1 first; missing axes start at 0.
         * @param size   Length of the window along each axis, NAXIS1 first; missing or 0 lengths, or lengths
         * going past the end of an axis, go to its end.
         * @return The new handle.
         * @throws FITSexception BAD_DIMEN if the window starts outside the view.
         */
        std::shared_ptr<ASDFarray> Window(const std::vector<size_t>& origin, const std::vector<size_t>& size) const;

        /**
         * @brief Make a handle showing one layer: one position along the last axis.
         * @details The new handle has one axis less. Nothing is read, as for Window().
         * @param index Position along the last axis, from 0.
         * @return The new handle.
         * @throws FITSexception BAD_DIMEN if the view has fewer than two axes or no such layer.
         */
        std::shared_ptr<ASDFarray> Layer(size_t index) const;

        /**
         * @brief Make an independent handle on the same array.
         * @details Pixels in memory, if any, are copied, with the changes made to them; a handle not loaded
         * gives a copy that isn't loaded either.
         * @return The new handle.
         */
        std::shared_ptr<ASDFarray> Clone() const;

        /**
         * @brief Make a handle on an image in memory.
         * @param cube The image; not null.
         * @return The new handle.
         * @throws FITSexception SHARED_NULPTR if @p cube is null.
         */
        static std::shared_ptr<ASDFarray> FromCube(std::shared_ptr<FITScube> cube);

#pragma endregion
    private:
        ASDFarray();                    //!< Empty handle, filled by asdfdetail::ArrayAccess
        friend struct asdfdetail::ArrayAccess;

        std::shared_ptr<const asdfdetail::ArraySource> fsource;  //!< Where the array is in its file; null in memory
        std::vector<uint64_t> forigin;  //!< First element of the view in the file, ASDF order, all axes; empty for no view
        std::vector<uint64_t> fsize;    //!< Size of the view in the file, ASDF order, all axes; empty for no view
        size_t fdropped = 0;            //!< Leading (slowest) axes of the view left out by Layer(); their size is 1
        std::shared_ptr<FITScube> fcube;//!< Pixels in memory, if loaded or set
        bool fmemoryOnly = false;       //!< The pixels in memory can't be read again from the file
        mutable std::string fcompression;   //!< Compression, once known
        mutable bool fcompressionKnown = false; //!< fcompression has been read
    };
}

#endif
