//
//  fits2asdf.h
//  asdf2fits
//
//  Written by Claude (Anthropic's AI model) for William Gillard.
//  Copyright (c) 2026 William Gillard.
//
//  SPDX-License-Identifier: BSD-3-Clause
//  Licence text: see the LICENSE file at the root of this project.

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
     * and tables are structured ndarrays. Data are written as physical values, so BSCALE, BZERO,
     * TSCALn and TZEROn are not carried over. Masked pixels become NaN (real images) or a BLANK
     * value (integer images: the lowest value of a signed type, the highest of an unsigned one).
     *
     * @param path    File to create.
     * @param replace Overwrite an existing file; otherwise an existing file throws FITSexception.
     * @throws FITSexception Also if a table has variable-length cells, which have no ASDF form.
     */
    void fits2asdf(const FITScube&,  const std::string& path, bool replace = false); //!< As the primary HDU
    void fits2asdf(const FITStable&, const std::string& path, bool replace = false); //!< After an empty primary HDU
    void fits2asdf(FITSmanager&,     const std::string& path, bool replace = false); //!< All HDUs of an open file
}

#endif
