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
     * @brief Load an image from an ndarray.
     * @details Booleans load as uint8_t and 16-bit reals as float. Pixels that are NaN, or equal to
     * BLANK in an integer image, are masked. The WCS is loaded if the header has CTYPE1, CRPIX1 or WCSAXES.
     * @param ndarray Ndarray giving the pixels, their type and the shape.
     * @param header  Header giving the keywords of the image.
     * @return The new image.
     * @throws FITSexception BAD_DATATYPE if the type has no image form, READ_ERROR if reading fails.
     */
    std::shared_ptr<FITScube> readImage(asdf_ndarray_t* ndarray, const FITShdu& header);

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
}

#endif
