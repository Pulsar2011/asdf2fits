//
//  fits2asdf.h
//  asdf2fits
//
//  Written by Claude (Anthropic's AI model) for William Gillard.
//  Copyright (c) 2026 William Gillard.
//
//  SPDX-License-Identifier: BSD-3-Clause
//  Licence text: see the LICENSE file at the root of this project.

/**
 * @file fits2asdf.h
 * @brief Writing of FITS content as ASDF files.
 */
#ifndef _DSL_fits2asdf_
#define _DSL_fits2asdf_

#include <string>

#include <DSTfits/FITSmanager.h>

namespace DSL
{
    /**
     * @brief Save FITS content as an ASDF file in the astropy FITS layout.
     * @details The file holds a value tagged tag:astropy.org:astropy/fits/fits-1.2.0, a list of
     * {header, data} HDUs that asdf-astropy reads back as an HDUList. Images are numeric ndarrays
     * and tables are structured ndarrays. Data are written as physical values, so the scaling of
     * the source isn't carried over: TSCALn and TZEROn are dropped, and BSCALE and BZERO are the
     * values FITSimg uses for the stored pixel type. Masked pixels become NaN (real images) or a BLANK
     * value (integer images: the lowest value of a signed type, the highest of an unsigned one).
     * CHECKSUM and DATASUM are dropped too.
     *
     * The image is written as the primary HDU.
     *
     * An image read from an ASDF file (an ASDFcube) is written as ASDFcube::WriteASDF() does instead: in
     * the layout of its file, with its whole tree and its companion arrays.
     *
     * @param cube    Image to save.
     * @param path    File to create.
     * @param replace Overwrite an existing file; otherwise an existing file throws FITSexception.
     * @throws FITSexception FILE_NOT_CREATED if the file exists and @p replace is false, or if it
     * can't be written. BAD_DATATYPE if the image holds no known type, MEMORY_ALLOCATION if a buffer
     * can't be allocated, WRITE_ERROR if the data can't be serialized.
     */
    void fits2asdf(const FITScube& cube, const std::string& path, bool replace = false);

    /**
     * @brief Save a table as an ASDF file in the astropy FITS layout, after an empty primary HDU.
     * @details The file has the layout described for
     * fits2asdf(const FITScube&, const std::string&, bool). Strings are as wide as the longest value.
     *
     * @param table   Table to save.
     * @param path    File to create.
     * @param replace Overwrite an existing file; otherwise an existing file throws FITSexception.
     * @throws FITSexception FILE_NOT_CREATED if the file exists and @p replace is false, or if it
     * can't be written. BAD_TFORM if a table has variable-length or empty cells, which have no ASDF
     * form, BAD_DIMEN if a column has the wrong number of rows, MEMORY_ALLOCATION if a buffer can't
     * be allocated, WRITE_ERROR if the data can't be serialized.
     */
    void fits2asdf(const FITStable& table, const std::string& path, bool replace = false);

    /**
     * @brief Save all HDUs of an open FITS file as an ASDF file in the astropy FITS layout.
     * @details The file has the layout described for
     * fits2asdf(const FITScube&, const std::string&, bool). An image HDU without dimensions is
     * written with its header only.
     * @note This moves the current HDU of @p fits.
     *
     * @param fits    Open FITS file.
     * @param path    File to create.
     * @param replace Overwrite an existing file; otherwise an existing file throws FITSexception.
     * @throws FITSexception FILE_NOT_CREATED if the file exists and @p replace is false, or if it
     * can't be written. Otherwise as the overloads for one image and one table, or anything that
     * reading from @p fits throws.
     */
    void fits2asdf(FITSmanager& fits, const std::string& path, bool replace = false);
}

#endif
