//
//  ASDFcube.h
//  asdf2fits
//
//  Written by Claude (Anthropic's AI model) for William Gillard.
//  Copyright (c) 2026 William Gillard.
//
//  SPDX-License-Identifier: BSD-3-Clause
//  Licence text: see the LICENSE file at the root of this project.

/**
 * @file ASDFcube.h
 * @brief Images that keep the ASDF file they come from: DSL::ASDFcube and DSL::ASDFimg\<T\>.
 */
#ifndef _DSL_ASDFcube_
#define _DSL_ASDFcube_

#include <cstdint>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <DSTfits/FITSimg.h>

#include "ASDFarray.h"
#include "ASDFgwcs.h"
#include "ASDFnode.h"

namespace DSL
{
    namespace asdfdetail
    {
        struct HeaderLink;
        struct CubeAccess;
        struct GwcsCache;
    }

    /** @brief Options of ASDFcube::WriteASDF(). */
    struct ASDFwriteOptions
    {
        bool replace = false;               //!< Overwrite an existing file; otherwise an existing file throws FITSexception
        /**
         * @brief Compression of the arrays written from memory: the image, and the arrays loaded or cut.
         * @details "keep" (the default) gives each array the compression it had in its file, "" none, or
         * "zlib", "bzp2" or "lz4". Arrays that weren't loaded are copied as they are, compressed or not.
         */
        std::string compression = "keep";
        std::string history;                //!< Description of a history entry added to the file; empty for none
    };

    /**
     * @class DSL::ASDFcube ASDFcube.h "asdf2fits/ASDFcube.h"
     * @brief The ASDF side of an image read from an ASDF file: the whole tree of the file, and the arrays
     * that go with the image.
     * @details An ASDFcube is never made alone: ASDFmanager::GetImage() returns an ASDFimg\<T\>, which is both
     * a FITSimg\<T\> and an ASDFcube. Code that knows only FITScube reaches this side with
     * `std::dynamic_pointer_cast<ASDFcube>(cube)`.
     *
     * **Tree and header.** Meta() is the tree of the whole file; the pixels of the image are at ImagePath().
     * The header of the image (FITScube::HDU()) is a flat view of the tree, which is what FITS tools read and
     * change. The tree is the reference for ASDF: WriteASDF() first copies into it what was changed in the
     * header since it was read (see SyncFromHDU()). Set() and Erase() change the tree and the header together;
     * after changing Meta() directly, call RefreshHDU().
     *
     * **Companions.** When the image is inside a mapping, as in the data models of Roman or JWST, the other
     * ndarrays of that mapping (err, dq, var_poisson...) are its companions. They are handles that read
     * nothing until asked (see ASDFarray). Window(), Layer() and Resize() give the companions of the same
     * shape as the image the same cut, as views that still read nothing. Rebin(), Overlay(), AddLayer() and
     * any other change of shape drop them, as their pixels can't be combined like the image's; a history
     * entry says so.
     *
     * **GWCS.** The first GWCS object found next to the image is its world coordinate system (see GWCS()):
     * ASDFimg\<T\> converts pixels to world coordinates and back with it, in place of the FITS WCS of the header.
     * A window of a 2D image (Window(), Resize()) shifts the GWCS to the window. Any other change of shape
     * removes the GWCS objects found next to the image, as they would no longer match its pixels; the FITS WCS
     * of the header is kept up to date by FITSimg.
     *
     * **Writing.** WriteASDF() writes the whole tree in its own layout: the image from memory, loaded arrays
     * from memory, arrays not loaded copied from their file as they are, still compressed, without being
     * read. Values asdf2fits doesn't interpret, such as GWCS objects, units and times, are written back
     * unchanged. The history entries of the file are kept, and one is added for each change recorded by
     * AddHistory().
     *
     * Failures throw FITSexception.
     */
    class ASDFcube
    {
    public:
        virtual ~ASDFcube();

        /** @brief Get the FITS side of the image. @return The image. */
        virtual FITScube&       Cube() = 0;
        virtual const FITScube& Cube() const = 0;   //!< @copydoc Cube()

#pragma region * Tree
        /** @brief Get the tree of the file. @return The root of the tree. */
        const ASDFnode& Meta() const { return ftree; }
        /** @brief Get the tree of the file, to change it; call RefreshHDU() afterwards. @return The root of the tree. */
        ASDFnode&       Meta()       { return ftree; }

        /** @brief Get where the pixels of the image are in the tree. @return The path of the image's ndarray. */
        const std::string& ImagePath() const { return fimagePath; }

        /**
         * @brief Set a value of the tree, and update the header.
         * @param path  Path of the value, from the root of the tree; missing mappings are added.
         * @param value The value.
         * @throws FITSexception As ASDFnode::SetPath().
         */
        void Set(const std::string& path, ASDFnode value);

        /**
         * @brief Remove a value of the tree, and update the header.
         * @param path Path of the value, from the root of the tree.
         * @return True if it was there.
         */
        bool Erase(const std::string& path);

#pragma endregion
#pragma region * Companions
        /**
         * @brief Get the names of the companion arrays.
         * @return Keys of the ndarrays next to the image in its mapping, in file order; empty if the image isn't
         * in a mapping of its own.
         */
        std::vector<std::string> CompanionNames() const;

        /**
         * @brief Get a companion array.
         * @details The handle reads nothing until it is loaded; changes made to it are written by WriteASDF().
         * @param name Its key, such as "err" or "dq".
         * @return The handle, or null if there is no such companion.
         */
        std::shared_ptr<ASDFarray> Companion(const std::string& name) const;

        /**
         * @brief Mask the pixels flagged in a data quality array.
         * @details The array is loaded if it isn't, and released afterwards so that it costs no memory.
         * @param bits Flags that mask a pixel; all by default.
         * @param name Key of the data quality companion.
         * @return The number of pixels flagged.
         * @throws FITSexception KEY_NO_EXIST if there is no such companion, BAD_DIMEN if its shape isn't the
         * image's, BAD_DATATYPE if it doesn't hold integers.
         */
        size_t MaskFromDQ(uint64_t bits = ~uint64_t(0), const std::string& name = "dq");

#pragma endregion
#pragma region * GWCS
        /**
         * @brief Get the world coordinate system of the image.
         * @details The first GWCS object found next to the image (below the mapping holding it, as roman/meta/wcs
         * for an image at roman/data), read on the first call and kept until the tree changes. A GWCS that can't
         * be used, because a kind of transform isn't known (see ASDFtransform::Register()) or because its pixel
         * axes aren't those of the image, gives null, with a warning on the first call.
         * @return The GWCS; null if there is none or it can't be used.
         */
        std::shared_ptr<const ASDFgwcs> GWCS() const;

        /**
         * @brief Get where the GWCS of the image is in the tree.
         * @return The path of the GWCS object found by GWCS(), even if it can't be used; empty if there is none.
         */
        std::string GWCSPath() const;

#pragma endregion
#pragma region * Header and tree
        /**
         * @brief Copy into the tree what was changed in the header.
         * @details A changed keyword updates the value it was read from; a removed one removes it; a new one
         * is added next to the image (a dotted name A.B as nested mappings a/b). In the astropy FITS layout
         * the cards of the header are updated instead. The structure keywords (BITPIX, NAXISn, BSCALE...) are
         * left out, as they follow the pixels.
         */
        void SyncFromHDU();

        /**
         * @brief Read the header again from the tree, after changing Meta().
         * @details Calls SyncFromHDU() first, so that changes made to the header aren't lost.
         */
        void RefreshHDU();

#pragma endregion
#pragma region * History and writing
        /**
         * @brief Record a change, written as a history entry by WriteASDF().
         * @param text Description of the change.
         */
        void AddHistory(const std::string& text);

        /** @brief Get the changes recorded for the history. @return Their descriptions, oldest first. */
        const std::vector<std::string>& PendingHistory() const { return fhistory; }

        /**
         * @brief Write the image and its tree as an ASDF file.
         * @details The file is written next to @p path under a temporary name, then renamed, so writing over
         * the file the image was read from is safe. Nothing of this object is changed.
         * @param path    File to create.
         * @param options Options; see ASDFwriteOptions.
         * @throws FITSexception FILE_NOT_CREATED if the file exists and @c options.replace is false, or if it
         * can't be written; BAD_OPTION for an unknown compression; READ_ERROR if an array can't be read.
         */
        void WriteASDF(const std::string& path, const ASDFwriteOptions& options = {}) const;

#pragma endregion
    protected:
        ASDFcube();                                     //!< Empty state, given by the library
        ASDFcube(const ASDFcube& other);                //!< Independent copy: the companions are copied too
        ASDFcube& operator=(const ASDFcube& other);     //!< Independent copy: the companions are copied too

        /**
         * @brief Update the tree after the image was cut to a window.
         * @param before Axis lengths before the cut, NAXIS1 first.
         * @param origin First pixel of the window, NAXIS1 first; missing axes start at 0.
         * @param what   Description of the change, for the history.
         */
        void AfterWindow(const std::vector<size_t>& before, const std::vector<size_t>& origin, const std::string& what);

        /**
         * @brief Update the tree after one layer of the image was taken.
         * @param before Axis lengths before, NAXIS1 first.
         * @param index  Position of the layer along the last axis.
         * @param what   Description of the change, for the history.
         */
        void AfterLayer(const std::vector<size_t>& before, size_t index, const std::string& what);

        /**
         * @brief Update the tree after the shape of the image changed otherwise.
         * @param before Axis lengths before, NAXIS1 first.
         * @param what   Description of the change, for the history.
         */
        void AfterReshape(const std::vector<size_t>& before, const std::string& what);

        /**
         * @brief Update the tree after the pixels were replaced by a result made from them.
         * @details Companions of the image's shape are dropped; GWCS objects are kept if the shape didn't change.
         * @param before Axis lengths of the image the result was made from, NAXIS1 first.
         * @param what   Description of the change, for the history.
         */
        void AfterProcessing(const std::vector<size_t>& before, const std::string& what);

    private:
        friend struct asdfdetail::CubeAccess;

        std::vector<size_t> axesNow() const;            //!< Axis lengths of the image now
        std::vector<std::pair<std::string, std::shared_ptr<ASDFarray>>> companions(const ASDFnode& tree) const; //!< Companions in a tree, path and handle
        void alignTo(const std::vector<size_t>& before);//!< Drop the companions left behind by a change of shape made without the hooks
        void forgetGWCS() const;                        //!< The tree changed: read the GWCS again when asked

        ASDFnode ftree;                                 //!< Tree of the file
        std::string fimagePath;                         //!< Path of the image's ndarray
        bool fcompanions = false;                       //!< The image's siblings are its companions
        std::unique_ptr<asdfdetail::HeaderLink> flink;  //!< How the header maps to the tree
        FITShdu fsnapshot;                              //!< The header when it last matched the tree
        std::vector<std::string> fhistory;              //!< Changes to record as history entries
        std::vector<size_t> faligned;                   //!< Axis lengths of the image the companions match
        std::unique_ptr<asdfdetail::GwcsCache> fgwcs;   //!< The GWCS, once read
    };

    /**
     * @class DSL::ASDFimg ASDFcube.h "asdf2fits/ASDFcube.h"
     * @brief An image read from an ASDF file: a FITSimg\<T\> that also keeps the file's tree and companions.
     * @details Everything that takes a FITScube or a FITSimg\<T\> takes it. Window(), Layer(), Rebin() and
     * Overlay() return ASDFimg\<T\> images, which keep the tree; see ASDFcube for what happens to the
     * companions and the GWCS. Results that other code builds as plain FITSimg\<T\> images get the tree back
     * with asdfLike().
     *
     * Instantiated for int8_t, uint8_t, int16_t, uint16_t, int32_t, uint32_t, int64_t, uint64_t, float and double.
     * @tparam T Pixel type.
     */
    template<typename T>
    class ASDFimg : public FITSimg<T>, public ASDFcube
    {
    public:
        /**
         * @brief Make an empty image, with no tree.
         * @param axes Axis lengths, NAXIS1 first.
         */
        explicit ASDFimg(const std::vector<size_t>& axes);
        ASDFimg(const ASDFimg<T>& other);                       //!< Copy: pixels, header, tree and companions
        ASDFimg<T>& operator=(const ASDFimg<T>& other);         //!< Copy: pixels, header, tree and companions
        ~ASDFimg() override;

        FITScube&       Cube() override       { return *this; } //!< @copydoc ASDFcube::Cube()
        const FITScube& Cube() const override { return *this; } //!< @copydoc ASDFcube::Cube()

        using FITScube::Rebin;
        using FITScube::Resize;
        using FITScube::WorldCoordinates;
        using FITScube::WorldCoordinatesVector;

        /**
         * @brief Get the world coordinates of a pixel.
         * @details With the GWCS of the image (see ASDFcube::GWCS()) when it has one and @p wcsIndex is 0;
         * otherwise as FITSimg, with the FITS WCS of the header. The other forms of WorldCoordinates() and
         * WorldCoordinatesVector() come here too.
         * @param pixel    Position, 0-based, NAXIS1 first.
         * @param wcsIndex Index of the WCS; 0 is the GWCS if there is one.
         * @return The world coordinates; NaN outside the bounding box of the GWCS.
         */
        worldCoords  WorldCoordinates(const pixelCoords& pixel, const int& wcsIndex = 0) const override;
        worldVectors WorldCoordinatesVector(const pixelVectors& pixels, const int& wcsIndex = 0) const override;   //!< As WorldCoordinates(), for several pixels
        /**
         * @brief Get the pixel of world coordinates.
         * @details As WorldCoordinates(): with the GWCS when there is one and @p wcsIndex is 0.
         * @param world    World coordinates.
         * @param wcsIndex Index of the WCS; 0 is the GWCS if there is one.
         * @return The position, 0-based, NAXIS1 first; NaN if it is outside the bounding box of the GWCS.
         */
        pixelCoords  World2Pixel(const worldCoords& world, const int& wcsIndex = 0) const override;
        pixelVectors World2PixelVector(const worldVectors& world, const int& wcsIndex = 0) const override;         //!< As World2Pixel(), for several positions

        /** @brief Tell whether the image has world coordinates. @return True with a GWCS, or as FITSimg with a FITS WCS. */
        bool HasWorldCoordinates() const override;
        /**
         * @brief Get the name of each world axis.
         * @param wcsIndex Index of the WCS; 0 is the GWCS if there is one.
         * @return The axes_names of the output frame of the GWCS, such as {"lon", "lat"}; otherwise as FITSimg (CTYPEi).
         */
        std::vector<std::string> WorldAxisNames(const int& wcsIndex = 0) const override;
        /**
         * @brief Get the unit of each world axis.
         * @param wcsIndex Index of the WCS; 0 is the GWCS if there is one.
         * @return The units of the output frame of the GWCS, such as {"deg", "deg"}; otherwise as FITSimg (CUNITi).
         */
        std::vector<std::string> WorldAxisUnits(const int& wcsIndex = 0) const override;

        std::shared_ptr<FITScube> Layer(const size_t& index) const override;                                   //!< As FITSimg; companions become layers too
        std::shared_ptr<FITScube> Window(size_t xMin, size_t yMin, size_t xSize, size_t ySize) const override;  //!< As FITSimg; companions become windows too
        std::shared_ptr<FITScube> Rebin(const std::vector<size_t>& nbin, bool doMean = false) const override;  //!< As FITSimg; companions of the image's shape are dropped
        std::shared_ptr<FITScube> Overlay(const FITScube::overlay& method = FITScube::overlay::mean,
                                          const std::pair<double,double>& clip = std::pair<double,double>(-1.,-1.)) const override;   //!< As FITSimg; companions of the image's shape are dropped
        void Resize(const size_t& xMin, const size_t& yMin, const size_t& xSize, const size_t& ySize) override; //!< As FITSimg; companions are cut too
        void AddLayer(const FITScube& layer) override;          //!< As FITSimg; companions of the image's shape are dropped
        void AddLayer(const FITSimg<T>& layer);                 //!< As FITSimg; companions of the image's shape are dropped

        /**
         * @brief Make an ASDF image from a plain one, with the tree of another.
         * @details The pixels are taken from @p image, which is left empty; its mask and header are copied.
         * @param image           Image giving the pixels, the mask and the header.
         * @param model           Image giving the tree, the companions and the history.
         * @param inheritKeywords Also copy the keywords of the model's header that @p image lacks, structure
         * keywords aside: for an image built with a header of its own, whose missing keywords weren't removed
         * on purpose and must not be removed from the tree.
         * @return The new image.
         */
        static std::shared_ptr<ASDFimg<T>> Adopt(FITSimg<T>& image, const ASDFcube& model, bool inheritKeywords = false);
    };

    extern template class ASDFimg<int8_t>;
    extern template class ASDFimg<uint8_t>;
    extern template class ASDFimg<int16_t>;
    extern template class ASDFimg<uint16_t>;
    extern template class ASDFimg<int32_t>;
    extern template class ASDFimg<uint32_t>;
    extern template class ASDFimg<int64_t>;
    extern template class ASDFimg<uint64_t>;
    extern template class ASDFimg<float>;
    extern template class ASDFimg<double>;

    /**
     * @brief Give a result built from an ASDF image the tree of that image.
     * @details For results that code knowing only FITScube builds as plain FITSimg\<T\> images, such as a
     * stack or a filtered image. Companions of the model's shape are dropped, as they don't describe the new
     * pixels; GWCS objects are kept if the shape is the model's. A history entry records @p what.
     * @param result The result. If nothing else refers to it, its pixels are moved rather than copied.
     * @param model  The image the result was made from.
     * @param what   Description of how the result was made, for the history.
     * @return An ASDFimg\<T\> holding the result; @p result itself if it is already an ASDF image.
     * @throws FITSexception BAD_DATATYPE if the result holds no known pixel type, SHARED_NULPTR if it is null.
     */
    std::shared_ptr<FITScube> asdfLike(std::shared_ptr<FITScube> result, const ASDFcube& model, const std::string& what = "processed");

    /**
     * @brief Write a tree as an ASDF file.
     * @details As ASDFcube::WriteASDF(), for a tree with no image of its own: one from ASDFmanager::GetTree()
     * after changing its metadata, or one built from scratch with ASDFarray::FromCube() handles. Arrays
     * loaded or made in memory are written from memory, the others are copied from their file as they are.
     * @param tree    Root of the tree.
     * @param path    File to create.
     * @param options Options; see ASDFwriteOptions.
     * @throws FITSexception As ASDFcube::WriteASDF().
     */
    void WriteASDF(const ASDFnode& tree, const std::string& path, const ASDFwriteOptions& options = {});
}

#endif
