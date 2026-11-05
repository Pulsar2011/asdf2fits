//
//  ASDFmanager.h
//  asdf2fits
//
//  Written by Claude (Anthropic's AI model) for William Gillard.
//  Copyright (c) 2026 William Gillard.
//
//  SPDX-License-Identifier: BSD-3-Clause
//  Licence text: see the LICENSE file at the root of this project.

#ifndef _DSL_ASDFmanager_
#define _DSL_ASDFmanager_

#include <memory>
#include <string>

#include <DSTfits/FITShdu.h>
#include <DSTfits/FITSimg.h>
#include <DSTfits/FITStable.h>

namespace DSL
{
    /**
     * @class DSL::ASDFmanager ASDFmanager.h "asdf2fits/ASDFmanager.h"
     * @brief Opens an ASDF file and hands out its content as DeepSkyFits objects.
     * @details Mirrors FITSmanager: HDU indices start at 1, names are matched without regard to
     * case, and failures throw FITSexception. Image and table data are only read from the file
     * when requested, so GetHeader() and friends stay cheap whatever the file size.
     *
     * Two kinds of file are understood:
     *  - the astropy FITS layout (a value tagged .../fits/fits-*), written by fits2asdf(): every
     *    list entry is one HDU, with its exact FITS header;
     *  - any other tree: each numeric ndarray is an image, each structured ndarray or
     *    .../table-* value is a table. Scalars are turned into keywords, nested keys being joined
     *    as META.INSTRUMENT.NAME. An HDU is named after its tree path, and the one called PRIMARY
     *    (or else the first found) is the primary HDU. GWCS objects are not converted.
     *
     * Values in an ASDF file are physical: BSCALE and BZERO found in a header are not applied.
     * NaN (real images) or the BLANK value (integer images) mark masked pixels.
     */
    class ASDFmanager
    {
    public:
        explicit ASDFmanager(const std::string& path);      //!< Open an existing ASDF file
        ~ASDFmanager();

        ASDFmanager(const ASDFmanager&) = delete;
        ASDFmanager& operator=(const ASDFmanager&) = delete;

        int NumberOfHeader() const;                         //!< Total number of HDU
        bool isOpen() const;
        void Close();                                       //!< Objects already handed out stay valid

        std::shared_ptr<FITShdu>   GetPrimaryHeader();
        std::shared_ptr<FITShdu>   GetHeaderAtIndex(const int&);
        std::shared_ptr<FITShdu>   GetHeader(const std::string&);

        std::shared_ptr<FITScube>  GetPrimary();
        std::shared_ptr<FITScube>  GetImageAtIndex(const int&);
        std::shared_ptr<FITScube>  GetImage(const std::string&);

        std::shared_ptr<FITStable> GetTableAtIndex(const int&);
        std::shared_ptr<FITStable> GetTable(const std::string&);

    private:
        struct Impl;
        std::unique_ptr<Impl> impl;
    };
}

#endif
