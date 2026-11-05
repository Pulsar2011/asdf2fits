//
//  ASDFmanager.h
//  asdf2fits
//
//  Written by Claude (Anthropic's AI model) for William Gillard.
//  Copyright (c) 2026 William Gillard.
//
//  SPDX-License-Identifier: BSD-3-Clause
//  Licence text: see the LICENSE file at the root of this project.

/**
 * @file ASDFmanager.h
 * @brief Reading of ASDF files as DeepSkyFits objects.
 */
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
     *
     * Each public call locks a mutex of the instance. Objects handed out are independent copies:
     * changing them doesn't affect the manager, and they outlive it.
     */
    class ASDFmanager
    {
    public:
        /**
         * @brief Open an existing ASDF file and list its HDUs.
         * @param path File to open.
         * @throws FITSexception FILE_NOT_OPENED if the file can't be opened or has no tree. BAD_DATATYPE,
         * BAD_DIMEN, BAD_TFORM or NOT_TABLE if the description of an HDU can't be read.
         */
        explicit ASDFmanager(const std::string& path);

        /** @brief Close the file. */
        ~ASDFmanager();

        ASDFmanager(const ASDFmanager&) = delete;               //!< Not copyable: the manager owns the file
        ASDFmanager& operator=(const ASDFmanager&) = delete;    //!< Not copyable: the manager owns the file

        /**
         * @brief Get the number of HDUs.
         * @return The total number of HDUs, 0 once the file is closed.
         */
        int NumberOfHeader() const;

        /**
         * @brief Tell whether the file is open.
         * @return False once Close() has been called.
         */
        bool isOpen() const;

        /**
         * @brief Close the file and forget its HDUs. Objects already handed out stay valid.
         * @details Calling it again has no effect. Afterwards, any call to get an HDU throws FITSexception.
         */
        void Close();

        /**
         * @brief Get a copy of the header of the primary HDU.
         * @return The header of HDU 1.
         * @throws FITSexception FILE_NOT_OPENED if the file is closed.
         */
        std::shared_ptr<FITShdu>   GetPrimaryHeader();

        /**
         * @brief Get a copy of the header of an HDU.
         * @param index Index of the HDU, starting at 1.
         * @return The header.
         * @throws FITSexception FILE_NOT_OPENED if the file is closed, BAD_HDU_NUM if there is no such HDU.
         */
        std::shared_ptr<FITShdu>   GetHeaderAtIndex(const int& index);

        /**
         * @brief Get a copy of the header of an HDU.
         * @param name Name of the HDU, matched without regard to case.
         * @return The header.
         * @throws FITSexception FILE_NOT_OPENED if the file is closed, BAD_HDU_NUM if there is no such HDU.
         */
        std::shared_ptr<FITShdu>   GetHeader(const std::string& name);

        /**
         * @brief Read the image of the primary HDU.
         * @return The image of HDU 1.
         * @throws FITSexception FILE_NOT_OPENED if the file is closed, NOT_IMAGE if HDU 1 is a table or
         * holds no data, BAD_DATATYPE or READ_ERROR if the pixels can't be loaded.
         */
        std::shared_ptr<FITScube>  GetPrimary();

        /**
         * @brief Read the image of an HDU.
         * @param index Index of the HDU, starting at 1.
         * @return The image.
         * @throws FITSexception FILE_NOT_OPENED if the file is closed, BAD_HDU_NUM if there is no such HDU,
         * NOT_IMAGE if it is a table or holds no data, BAD_DATATYPE or READ_ERROR if the pixels can't be loaded.
         */
        std::shared_ptr<FITScube>  GetImageAtIndex(const int& index);

        /**
         * @brief Read the image of an HDU.
         * @param name Name of the HDU, matched without regard to case.
         * @return The image.
         * @throws FITSexception FILE_NOT_OPENED if the file is closed, BAD_HDU_NUM if there is no such HDU,
         * NOT_IMAGE if it is a table or holds no data, BAD_DATATYPE or READ_ERROR if the pixels can't be loaded.
         */
        std::shared_ptr<FITScube>  GetImage(const std::string& name);

        /**
         * @brief Read the table of an HDU.
         * @param index Index of the HDU, starting at 1.
         * @return The table.
         * @throws FITSexception FILE_NOT_OPENED if the file is closed, BAD_HDU_NUM if there is no such HDU,
         * NOT_TABLE if it isn't a table, BAD_DIMEN, BAD_TFORM or READ_ERROR if its columns can't be loaded.
         */
        std::shared_ptr<FITStable> GetTableAtIndex(const int& index);

        /**
         * @brief Read the table of an HDU.
         * @param name Name of the HDU, matched without regard to case.
         * @return The table.
         * @throws FITSexception FILE_NOT_OPENED if the file is closed, BAD_HDU_NUM if there is no such HDU,
         * NOT_TABLE if it isn't a table, BAD_DIMEN, BAD_TFORM or READ_ERROR if its columns can't be loaded.
         */
        std::shared_ptr<FITStable> GetTable(const std::string& name);

    private:
        struct Impl;                //!< Open file and list of HDUs
        std::unique_ptr<Impl> impl; //!< Private state of the manager
    };
}

#endif
