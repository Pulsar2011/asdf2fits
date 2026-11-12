//
//  ASDFgwcs.h
//  asdf2fits
//
//  Written by Claude (Anthropic's AI model) for William Gillard.
//  Copyright (c) 2026 William Gillard.
//
//  SPDX-License-Identifier: BSD-3-Clause
//  Licence text: see the LICENSE file at the root of this project.

/**
 * @file ASDFgwcs.h
 * @brief World coordinates of an ASDF image: DSL::ASDFgwcs, a GWCS object read from the tree.
 */
#ifndef _DSL_ASDFgwcs_
#define _DSL_ASDFgwcs_

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "ASDFnode.h"
#include "ASDFtransform.h"

namespace DSL
{
    /** @brief One frame of a GWCS: where the coordinates are between two steps. */
    struct ASDFframe
    {
        std::string name;                   //!< Name, such as "detector", "v2v3" or "world"
        std::string kind;                   //!< Kind, from its tag: "gwcs/frame2d", "gwcs/celestial_frame"...; empty if untagged
        std::vector<std::string> axesNames; //!< Name of each axis; may be empty
        std::vector<std::string> units;     //!< Unit of each axis, such as "pixel", "arcsec" or "deg"; may be empty
        std::string referenceFrame;         //!< Celestial reference frame, such as "icrs"; empty if none

        /** @brief Tell whether the frame is celestial. @return True if its first two axes are a longitude and a latitude. */
        bool IsCelestial() const { return kind == "gwcs/celestial_frame"; }
    };

    /**
     * @class DSL::ASDFgwcs ASDFgwcs.h "asdf2fits/ASDFgwcs.h"
     * @brief A GWCS: a chain of frames, from the pixels to the sky, linked by transforms.
     * @details Read from a node tagged "tag:stsci.edu:gwcs/wcs-...", such as roman/meta/wcs in a Roman file. The
     * steps of the node give the frames, from the first (the pixels) to the last (the world), and the transform
     * from each frame to the next. Transforms are read by ASDFtransform::Read(): a kind of transform the library
     * doesn't know makes Read() fail, until a builder is registered for it.
     *
     * Pixel coordinates start at 0 at the centre of the first pixel, as in GWCS (and FITScube): a FITS pixel
     * position is one more. The bounding box, when the file gives one, is where the transforms are valid: as in
     * GWCS, a pixel outside gives NaN, and so does a world position whose pixel falls outside.
     *
     * World to pixel uses the inverses of the transforms: those given by the file, or derived from the
     * transforms themselves (see ASDFtransform::Inverse()). If a step has neither, the pixel is found by
     * Newton's method on the transforms.
     *
     * An ASDFgwcs is immutable once read, so it can be shared and used from several threads.
     * Failures throw WCSexception.
     */
    class ASDFgwcs
    {
    public:
        /**
         * @brief Read a GWCS.
         * @param wcs Node of the GWCS object: a mapping with "steps".
         * @return The GWCS.
         * @throws WCSexception BAD_OPTION if the node isn't a GWCS, or if a transform can't be read; the message
         * names the kind of transform that isn't known.
         */
        static std::shared_ptr<const ASDFgwcs> Read(const ASDFnode& wcs);

        /** @brief Get the name of the GWCS. @return The name given in the file, such as "FIT-LVL2-GAIADR3"; may be empty. */
        const std::string& Name() const { return fname; }

        /** @brief Get the frames. @return The frames, from the pixels to the world. */
        const std::vector<ASDFframe>& Frames() const { return fframes; }

        /** @brief Get the first frame. @return The frame of the pixels. */
        const ASDFframe& InputFrame() const { return fframes.front(); }

        /** @brief Get the last frame. @return The frame of the world coordinates. */
        const ASDFframe& OutputFrame() const { return fframes.back(); }

        /** @brief Get the number of pixel axes. @return The number of inputs of the first transform. */
        size_t PixelAxes() const;

        /** @brief Get the number of world axes. @return The number of outputs of the last transform. */
        size_t WorldAxes() const;

        /**
         * @brief Get the transform from one frame to the next.
         * @param step Index of the first frame, from 0.
         * @return The transform from frame @p step to frame @p step + 1.
         * @throws WCSexception BAD_DIMEN if there is no such step.
         */
        const std::shared_ptr<const ASDFtransform>& Step(size_t step) const;

        /**
         * @brief Get the bounding box.
         * @return The range [lower, upper] of each pixel axis, NAXIS1 first; empty if the file gives none.
         */
        const std::vector<std::pair<double,double>>& BoundingBox() const { return fbox; }

        /** @brief Tell whether world to pixel is computed by inverse transforms. @return False if it is found numerically. */
        bool HasInverse() const;

        /**
         * @brief Convert a pixel position to world coordinates.
         * @param pixel           Position, 0-based, NAXIS1 first.
         * @param withBoundingBox Give NaN outside the bounding box.
         * @return The world coordinates, in the units of the output frame (degrees on the sky).
         * @throws WCSexception BAD_DIMEN if @p pixel hasn't PixelAxes() values.
         */
        std::vector<double> PixelToWorld(const std::vector<double>& pixel, bool withBoundingBox = true) const;

        /**
         * @brief Convert world coordinates to a pixel position.
         * @param world           Coordinates, in the units of the output frame.
         * @param withBoundingBox Give NaN when the pixel is outside the bounding box.
         * @return The position, 0-based, NAXIS1 first; NaN if it can't be found.
         * @throws WCSexception BAD_DIMEN if @p world hasn't WorldAxes() values.
         */
        std::vector<double> WorldToPixel(const std::vector<double>& world, bool withBoundingBox = true) const;

        /**
         * @brief Convert coordinates from one frame to another, forward or backward.
         * @details No bounding box is applied.
         * @param from   Name of the frame of @p coords.
         * @param to     Name of the frame wanted.
         * @param coords Coordinates in frame @p from.
         * @return The coordinates in frame @p to; NaN if they can't be found.
         * @throws WCSexception BAD_OPTION if a frame isn't known, BAD_DIMEN if @p coords hasn't the right size.
         */
        std::vector<double> Convert(const std::string& from, const std::string& to, const std::vector<double>& coords) const;

    private:
        ASDFgwcs() = default;   //!< Filled by Read()

        size_t frameIndex(const std::string& name) const;                                      //!< Index of a frame by name
        std::vector<double> run(size_t from, size_t to, std::vector<double> coords) const;      //!< Between two frames by index
        std::vector<double> solve(size_t from, size_t to, const std::vector<double>& target) const; //!< Backward, by Newton's method
        bool inside(const std::vector<double>& pixel) const;                                     //!< Inside the bounding box

        std::string fname;                                          //!< Name of the GWCS
        std::vector<ASDFframe> fframes;                             //!< Frames, pixels first
        std::vector<std::shared_ptr<const ASDFtransform>> fsteps;   //!< Transform from frame k to frame k+1
        std::vector<std::shared_ptr<const ASDFtransform>> finverses;//!< Inverse of each step; null if none
        std::vector<std::pair<double,double>> fbox;                 //!< Bounding box of the pixels; empty if none
    };
}

#endif
