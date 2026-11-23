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
#include <vector>

#include <DSTfits/FITShdu.h>
#include <DSTfits/FITSimg.h>
#include <DSTfits/FITStable.h>

#include "ASDFarray.h"
#include "ASDFcube.h"
#include "ASDFnode.h"

namespace DSL
{
    /** @brief Options of ASDFmanager. */
    struct ASDFreadOptions
    {
        /**
         * @brief Size in bytes above which compressed arrays are decompressed to a temporary file rather than
         * to memory; 0 for libasdf's default.
         */
        size_t maxMemoryBytes = 0;
        std::string tmpDir;     //!< Directory of those temporary files; empty for the system's
    };

    /** @brief Where an HDU of an ASDF file is in the tree of the file; see ASDFmanager::GetPathAtIndex(). */
    struct ASDFhduPath
    {
        /**
         * @brief Path of the subtree the HDU is read from.
         * @details In the astropy FITS layout, its entry of the list, as "fits/0"; in any other tree, the mapping
         * whose values are its keywords, as "roman". Empty when no subtree is its own: data directly under the
         * root, or the PRIMARY HDU made for a tree without data.
         */
        std::string node;
        std::string data;   //!< Path of its data value, as "roman/data" or "fits/0/data"; empty for the PRIMARY HDU made for a tree without data
    };

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
     *    (or else the first found) is the primary HDU. GWCS objects are not turned into keywords: images
     *    use them for their world coordinates (see ASDFcube::GWCS()).
     *
     * Values in an ASDF file are physical: BSCALE and BZERO found in a header are not applied.
     * NaN (real images) or the BLANK value (integer images) mark masked pixels.
     *
     * Images are ASDFimg\<T\>: FITSimg\<T\> images that also keep the tree of the file and the arrays
     * that go with them (see ASDFcube). Only their own pixels are read; the other arrays are read when asked.
     *
     * Each public call locks a mutex of the instance. Objects handed out are independent copies:
     * changing them doesn't affect the manager, and they outlive it. The arrays they haven't read yet
     * keep the file open until they are destroyed.
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

        /**
         * @brief Open an existing ASDF file and list its HDUs, with options.
         * @param path    File to open.
         * @param options How large compressed arrays are decompressed.
         * @throws FITSexception As ASDFmanager(const std::string&).
         */
        ASDFmanager(const std::string& path, const ASDFreadOptions& options);

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
         * Arrays handed out and not read yet keep the file open until they are destroyed.
         */
        void Close();

        /**
         * @brief Get the tree of the file.
         * @details No array is read: ndarrays are ASDFarray handles.
         * @return A copy of the whole tree.
         * @throws FITSexception FILE_NOT_OPENED if the file is closed.
         */
        ASDFnode GetTree();

        /**
         * @brief Get an ndarray of the file, without reading it.
         * @param path Path of the ndarray in the tree, as "roman/dq".
         * @return A handle that reads the array when asked.
         * @throws FITSexception FILE_NOT_OPENED if the file is closed, KEY_NO_EXIST if there is nothing at
         * @p path, NOT_IMAGE if it isn't an ndarray.
         */
        std::shared_ptr<ASDFarray> GetArray(const std::string& path);

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
         * @brief Get where an HDU is in the tree of the file (see GetTree()).
         * @details The paths are those of ASDFnode, so that GetTree().FindPath() finds the nodes they name.
         * @param index Index of the HDU, starting at 1.
         * @return The paths of the subtree the HDU is read from and of its data.
         * @throws FITSexception FILE_NOT_OPENED if the file is closed, BAD_HDU_NUM if there is no such HDU.
         */
        ASDFhduPath                GetPathAtIndex(const int& index);

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
         * @brief Read part of the image of an HDU.
         * @details Only the pixels of the window are read. The companions of the same shape show the same
         * window, a GWCS of a 2D image is shifted to the window (other GWCS objects next to the image are
         * removed), and a history entry records the cut (see ASDFcube). A FITS WCS of the header is shifted to
         * the window.
         * @param name   Name of the HDU, matched without regard to case.
         * @param origin First pixel of the window, NAXIS1 first; missing axes start at 0.
         * @param size   Length of the window along each axis, NAXIS1 first; missing or 0 lengths, or lengths
         * going past the end of an axis, go to its end.
         * @return The image of the window.
         * @throws FITSexception As GetImage(), or BAD_DIMEN if the window starts outside the image.
         */
        std::shared_ptr<FITScube>  GetImageWindow(const std::string& name, const std::vector<size_t>& origin, const std::vector<size_t>& size = {});

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
