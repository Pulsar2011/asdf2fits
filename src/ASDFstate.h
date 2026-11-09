//
//  ASDFstate.h
//  asdf2fits
//
//  Private: giving an image read by ASDFmanager its ASDF state.
//
//  Written by Claude (Anthropic's AI model) for William Gillard.
//  Copyright (c) 2026 William Gillard.
//
//  SPDX-License-Identifier: BSD-3-Clause
//  Licence text: see the LICENSE file at the root of this project.

/**
 * @file ASDFstate.h
 * @brief Access to the private state of ASDFcube, for ASDFmanager.
 * @note Private header: it isn't installed with the library.
 */
#ifndef _DSL_ASDFstate_
#define _DSL_ASDFstate_

#include <string>
#include <vector>

#include "ASDFcube.h"
#include "ASDFheader.h"

namespace DSL::asdfdetail
{
    /** @brief The ASDF state of an image as read from its file. */
    struct ImageState
    {
        ASDFnode tree;          //!< Tree of the file, independent of the manager's
        std::string imagePath;  //!< Path of the image's ndarray
        bool companions = false;//!< The image's siblings are its companions
        HeaderLink link;        //!< How the header maps to the tree; @c keyPaths is filled by attach()
    };

    /** @brief Access to the private state of ASDFcube. */
    struct CubeAccess
    {
        /**
         * @brief Give an image its ASDF state.
         * @details The header snapshot and the shape the companions match are taken from the image as it is.
         * @param cube  The image, as read.
         * @param state Its state.
         */
        static void attach(ASDFcube& cube, ImageState state);

        /**
         * @brief Update the tree of an image read as a window of its array.
         * @param cube   The image.
         * @param before Axis lengths of the whole array, NAXIS1 first.
         * @param origin First pixel of the window, NAXIS1 first.
         * @param what   Description of the change, for the history.
         */
        static void afterWindow(ASDFcube& cube, const std::vector<size_t>& before, const std::vector<size_t>& origin, const std::string& what)
        {
            // The companions still match the whole array
            cube.faligned = before;
            cube.AfterWindow(before, origin, what);
        }

        /**
         * @brief Update the tree of an image whose pixels were replaced by a result made from them.
         * @param cube   The image.
         * @param before Axis lengths of the image the result was made from, NAXIS1 first.
         * @param what   Description of the change, for the history.
         */
        static void afterProcessing(ASDFcube& cube, const std::vector<size_t>& before, const std::string& what)
        {
            cube.AfterProcessing(before, what);
        }
    };
}

#endif
