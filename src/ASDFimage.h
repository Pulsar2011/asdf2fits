//
//  ASDFimage.h
//  asdf2fits
//
//  Private: conversion between ndarrays and FITSimg<T>.
//
//  Licensed under PolyForm Noncommercial 1.0.0
//  You may use, share and adapt this code for noncommercial purposes.
//  Commercial use requires a separate licence from the author.
//  Licence text: https://polyformproject.org/licenses/noncommercial/1.0.0/

#ifndef _DSL_ASDFimage_
#define _DSL_ASDFimage_

#include <memory>

#include <DSTfits/FITSimg.h>

#include "ASDFhandle.h"

namespace DSL::asdfdetail
{
    /** @brief FITS axis lengths of an ndarray: NAXIS1 is its last, fastest varying, dimension. */
    std::vector<uint64_t> fitsAxes(const asdf_ndarray_t* ndarray);

    /** @brief True if images of this ndarray type can be loaded. */
    bool isImageType(asdf_scalar_datatype_t type);

    /** @brief Load an image; \c header provides its keywords, the ndarray its type and shape. */
    std::shared_ptr<FITScube> readImage(asdf_ndarray_t* ndarray, const FITShdu& header);

    /**
     * @brief Copy a cube into a new ndarray value of the file.
     * @param header Header of the output HDU; receives the structure keywords and BLANK if needed.
     */
    asdf_value_t* writeImage(asdf_file_t* file, const FITScube& cube, FITShdu& header);
}

#endif
