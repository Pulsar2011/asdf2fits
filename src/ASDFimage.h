//
//  ASDFimage.h
//  asdf2fits
//
//  Private: conversion between ndarrays and FITSimg<T>.
//
//  Written by Claude (Anthropic's AI model) for William Gillard.
//  Copyright (c) 2026 William Gillard.
//
//  SPDX-License-Identifier: BSD-3-Clause
//  Licence text: see the LICENSE file at the root of this project.

/**
 * @file ASDFimage.h
 * @brief Conversion between ndarrays and FITSimg\<T\>.
 * @note Private header: it isn't installed with the library.
 */
#ifndef _DSL_ASDFimage_
#define _DSL_ASDFimage_

#include <memory>

#include <DSTfits/FITSimg.h>

#include "ASDFhandle.h"
#include "ASDFnode.h"
#include "ASDFsource.h"

namespace DSL::asdfdetail
{
    /**
     * @brief Get the FITS axis lengths of an ndarray.
     * @param ndarray Ndarray to describe.
     * @return Axis lengths from NAXIS1 up: NAXIS1 is the last, fastest varying, ndarray dimension.
     */
    std::vector<uint64_t> fitsAxes(const asdf_ndarray_t* ndarray);

    /**
     * @brief Tell whether images of this ndarray type can be loaded.
     * @param type Element type of the ndarray.
     * @return True for integers of 8 to 64 bits, 16, 32 and 64-bit reals, and booleans.
     */
    bool isImageType(asdf_scalar_datatype_t type);

    /**
     * @brief Load an image, or part of one, from an ndarray.
     * @details Booleans load as uint8_t and 16-bit reals as float. Pixels that are NaN, or equal to
     * BLANK in an integer image, are masked. The WCS is loaded if the header has CTYPE1, CRPIX1 or WCSAXES.
     * @param ndarray   Ndarray giving the pixels, their type and the shape.
     * @param header    Header giving the keywords of the image.
     * @param tile      Part to read; the whole array by default.
     * @param asdfImage Make an ASDFimg\<T\> rather than a FITSimg\<T\>; its ASDF state is left empty.
     * @return The new image.
     * @throws FITSexception BAD_DATATYPE if the type has no image form, BAD_DIMEN if the tile isn't inside
     * the array, READ_ERROR if reading fails.
     */
    std::shared_ptr<FITScube> readImage(asdf_ndarray_t* ndarray, const FITShdu& header, const Tile& tile = {}, bool asdfImage = false);

    /**
     * @brief Copy a cube into a new ndarray value of the file.
     * @details Masked pixels are written as NaN in a real image. In an integer image they take the
     * lowest value of a signed type or the highest of an unsigned one, which is stored as BLANK.
     * @param file   File the value is created in.
     * @param cube   Cube to copy.
     * @param header Header of the output HDU; receives the structure keywords and BLANK if needed.
     * @return The new value.
     * @throws FITSexception BAD_DATATYPE if the cube holds no known type, MEMORY_ALLOCATION if the
     * pixel buffer can't be allocated, WRITE_ERROR if the value can't be created.
     */
    asdf_value_t* writeImage(asdf_file_t* file, const FITScube& cube, FITShdu& header);

    /** @brief Pixels of a cube, ready to be written as a block. */
    struct PixelBlock
    {
        const void* data = nullptr;         //!< The pixels; borrowed from the cube, or from @c owned
        size_t bytes = 0;                   //!< Size of the pixels
        std::vector<unsigned char> owned;   //!< Copy of the pixels, made when masked pixels had to be replaced
        asdf_scalar_datatype_t type = ASDF_DATATYPE_UNKNOWN;    //!< Element type
        std::vector<uint64_t> shape;        //!< Shape, slowest axis first
        ASDFnode blank;                     //!< Value given to masked pixels of an integer image; null if none
    };

    /**
     * @brief Get the pixels of a cube to write them as a block.
     * @details The pixels are borrowed from the cube, which must outlive the result, unless some are masked:
     * they are then copied, with NaN in a real image, and in an integer image the lowest value of a signed
     * type or the highest of an unsigned one.
     * @param cube Cube to write.
     * @return The pixels.
     * @throws FITSexception BAD_DATATYPE if the cube holds no known type.
     */
    PixelBlock cubePixels(const FITScube& cube);

    /**
     * @brief Get the YAML name of an element type.
     * @param type A type with an image form.
     * @return "float32", "uint16"...
     */
    std::string datatypeName(asdf_scalar_datatype_t type);

    /**
     * @brief Get the element type of a cube's pixels, without copying them.
     * @param cube The cube.
     * @return Its type, or ASDF_DATATYPE_UNKNOWN if it holds no known type.
     */
    asdf_scalar_datatype_t cubeType(const FITScube& cube);

    /**
     * @brief Copy a cube.
     * @param cube The cube; a FITSimg\<T\> of a known type.
     * @return A FITSimg\<T\> holding a copy of its pixels, mask and header.
     * @throws FITSexception BAD_DATATYPE if the cube holds no known type.
     */
    std::shared_ptr<FITScube> copyCube(const FITScube& cube);

    /**
     * @brief Copy part of a cube, with its mask.
     * @param cube    The cube.
     * @param origin  First pixel, one position per axis, NAXIS1 first.
     * @param size    Length along each axis, NAXIS1 first; the part must lie inside the cube.
     * @param dropped Last axes (NAXISn down) left out of the copy; their length must be 1.
     * @return A FITSimg\<T\> of the same type, named as the cube.
     * @throws FITSexception BAD_DATATYPE if the cube holds no known type.
     */
    std::shared_ptr<FITScube> cubeWindow(const FITScube& cube, const std::vector<size_t>& origin, const std::vector<size_t>& size, size_t dropped = 0);
}

#endif
