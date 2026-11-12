//
//  ASDFcube.cxx
//  asdf2fits
//
//  Written by Claude (Anthropic's AI model) for William Gillard.
//  Copyright (c) 2026 William Gillard.
//
//  SPDX-License-Identifier: BSD-3-Clause
//  Licence text: see the LICENSE file at the root of this project.

/**
 * @file ASDFcube.cxx
 * @brief Implementation of DSL::ASDFcube, DSL::ASDFimg\<T\> and DSL::asdfLike().
 */
#include "ASDFcube.h"

#include <algorithm>
#include <cstring>
#include <iostream>
#include <mutex>
#include <type_traits>

#include "ASDFheader.h"
#include "ASDFimage.h"
#include "ASDFstate.h"
#include "ASDFwrite.h"

namespace DSL
{
    using namespace asdfdetail;

    namespace
    {
        const char* const kWho = "ASDFcube";    //!< Name given in the exceptions thrown here

        /**
         * @brief Get the axis lengths of a cube.
         * @param cube The cube.
         * @return Its lengths, NAXIS1 first.
         */
        std::vector<size_t> axesOf(const FITScube& cube)
        {
            std::vector<size_t> axes;
            for(size_t i = 1; i <= cube.GetDimension(); i++)
                axes.push_back(cube.Size(i));
            return axes;
        }

        /**
         * @brief Write axis lengths for a history entry.
         * @param axes Lengths, NAXIS1 first.
         * @return The lengths joined by 'x', as "4088x4088".
         */
        std::string axesText(const std::vector<size_t>& axes)
        {
            std::string text;
            for(size_t length : axes)
                text += (text.empty() ? "" : "x") + std::to_string(length);
            return text.empty() ? std::string("0") : text;
        }

        /**
         * @brief Join names for a history entry.
         * @param names The names.
         * @return "a, b, c".
         */
        std::string join(const std::vector<std::string>& names)
        {
            std::string text;
            for(const std::string& name : names)
                text += (text.empty() ? "" : ", ") + name;
            return text;
        }

        /**
         * @brief Find the GWCS objects below a node: nodes tagged gwcs/..., not looking inside them.
         * @param top The node.
         * @return Their paths below @p top, in file order.
         */
        std::vector<std::string> gwcsObjects(const ASDFnode& top)
        {
            std::vector<std::string> found;
            top.Visit([&](const std::string& path, const ASDFnode& node)
            {
                if(node.HasTag("gwcs/") && !path.empty())
                {
                    found.push_back(path);
                    return false;
                }
                return true;
            });
            return found;
        }

        /**
         * @brief Find the tag prefix and version of the compositions, concatenations and shifts of a GWCS.
         * @details They are taken from the first of them in the GWCS. One version of the transform schemas gives
         * these three kinds the same version.
         * @param gwcs    The GWCS object.
         * @param prefix  Receives the tag up to the kind, as "tag:stsci.edu:asdf/transform/".
         * @param version Receives the version, as "1.3.0".
         * @return False if the GWCS has none of them.
         */
        bool transformTags(const ASDFnode& gwcs, std::string& prefix, std::string& version)
        {
            bool found = false;
            gwcs.Visit([&](const std::string&, const ASDFnode& node)
            {
                if(found)
                    return false;
                const std::string kind = ASDFtransform::KindOf(node.Tag());
                for(const char* name : {"compose", "concatenate", "shift"})
                    if(kind == std::string("transform/") + name)
                    {
                        const std::string& tag = node.Tag();
                        const size_t at = tag.rfind(std::string("/") + name + "-");
                        if(at != std::string::npos)
                        {
                            prefix = tag.substr(0, at + 1);
                            version = tag.substr(at + 2 + std::strlen(name));
                            found = true;
                        }
                        return false;
                    }
                return true;
            });
            return found;
        }

        /**
         * @brief Make a list of texts.
         * @param items The texts.
         * @return A sequence node.
         */
        ASDFnode textList(std::initializer_list<const char*> items)
        {
            ASDFnode list = ASDFnode::Sequence();
            for(const char* item : items)
                list.Append(ASDFnode(item));
            return list;
        }

        /**
         * @brief Make the GWCS of an image work on a window of it.
         * @details The first transform becomes a composition: a shift by the origin of the window on each axis,
         * then the transform as it was. The bounding box moves to the composition, shifted and cut to the
         * window, in the same form as it had (none becomes the whole window). The tags take the version of
         * the compositions, concatenations or shifts already in the GWCS, so that gwcs reads the result.
         * @param wcs    The GWCS object.
         * @param origin First pixel of the window, NAXIS1 first.
         * @param size   Lengths of the window, NAXIS1 first.
         * @return False if the GWCS can't be read, isn't two-dimensional, or has no tag to copy; it is then
         * left as it was.
         */
        bool shiftGWCS(ASDFnode& wcs, const std::vector<size_t>& origin, const std::vector<size_t>& size)
        {
            if(origin.size() != 2 || size.size() != 2)
                return false;
            std::shared_ptr<const ASDFgwcs> parsed;
            try
            {
                parsed = ASDFgwcs::Read(wcs);
            }
            catch(const FITSexception&)
            {
                return false;
            }
            std::string prefix, version;
            if(parsed->PixelAxes() != 2 || !transformTags(wcs, prefix, version))
                return false;

            ASDFnode& first = wcs.At("steps/0");
            ASDFnode old = first.At("transform");
            ASDFnode oldBox;
            if(const ASDFnode* box = old.Find("bounding_box"))
            {
                oldBox = *box;
                old.Erase("bounding_box");
            }

            // The window's pixel (i, j) is the image's pixel (i + x0, j + y0)
            ASDFnode shifts = ASDFnode::Mapping(prefix + "concatenate-" + version);
            ASDFnode& parts = shifts["forward"];
            for(size_t k = 0; k < 2; k++)
            {
                ASDFnode shift = ASDFnode::Mapping(prefix + "shift-" + version);
                shift["inputs"] = textList({"x"});
                shift["offset"] = static_cast<double>(origin[k]);
                shift["outputs"] = textList({"y"});
                parts.Append(std::move(shift));
            }
            shifts["inputs"] = textList({"x0", "x1"});
            shifts["outputs"] = textList({"y0", "y1"});

            std::vector<std::pair<double,double>> ranges;
            const std::vector<std::pair<double,double>>& before = parsed->BoundingBox();
            for(size_t k = 0; k < 2; k++)
            {
                const double lower = -0.5, upper = static_cast<double>(size[k]) - 0.5;
                if(before.empty())
                    ranges.emplace_back(lower, upper);
                else
                    ranges.emplace_back(std::max(before[k].first - static_cast<double>(origin[k]), lower),
                                        std::min(before[k].second - static_cast<double>(origin[k]), upper));
            }

            ASDFnode compose = ASDFnode::Mapping(prefix + "compose-" + version);
            if(oldBox.IsMapping())
            {
                ASDFnode box = ASDFnode::Mapping(oldBox.Tag());
                box["ignore"] = ASDFnode::Sequence();
                ASDFnode& intervals = box["intervals"];
                for(size_t k = 0; k < 2; k++)
                {
                    ASDFnode range = ASDFnode::Sequence();
                    range.Append(ranges[k].first);
                    range.Append(ranges[k].second);
                    intervals["x" + std::to_string(k)] = std::move(range);
                }
                box["order"] = ASDFnode("F");
                compose["bounding_box"] = std::move(box);
            }
            else
            {
                // The older form: a list of [lower, upper], the last axis first
                ASDFnode box = ASDFnode::Sequence();
                for(size_t k = 2; k-- > 0;)
                {
                    ASDFnode range = ASDFnode::Sequence();
                    range.Append(ranges[k].first);
                    range.Append(ranges[k].second);
                    box.Append(std::move(range));
                }
                compose["bounding_box"] = std::move(box);
            }

            ASDFnode& chain = compose["forward"];
            chain.Append(std::move(shifts));
            chain.Append(old);
            compose["inputs"] = textList({"x0", "x1"});
            if(const ASDFnode* outputs = old.Find("outputs"))
                compose["outputs"] = *outputs;
            first["transform"] = std::move(compose);

            if(ASDFnode* shape = wcs.Find("pixel_shape"); shape && !shape->IsNull())
            {
                ASDFnode lengths = ASDFnode::Sequence();
                lengths.Append(static_cast<int64_t>(size[0]));
                lengths.Append(static_cast<int64_t>(size[1]));
                *shape = std::move(lengths);
            }
            return true;
        }

        /**
         * @brief Shift the GWCS objects found below a node to a window, or remove those that can't be.
         * @param tree    Root of the tree.
         * @param scope   Path of the node to search; empty for the whole tree.
         * @param origin  First pixel of the window, NAXIS1 first; empty to remove them all.
         * @param size    Lengths of the window, NAXIS1 first.
         * @param shifted Receives the paths of the objects shifted.
         * @return Paths of the objects removed.
         */
        std::vector<std::string> fitGWCS(ASDFnode& tree, const std::string& scope, const std::vector<size_t>& origin,
                                         const std::vector<size_t>& size, std::vector<std::string>& shifted)
        {
            ASDFnode* top = tree.FindPath(scope);
            if(!top)
                return {};

            std::vector<std::string> doomed;
            for(const std::string& path : gwcsObjects(*top))
            {
                ASDFnode& wcs = top->At(path);
                if(!origin.empty() && ASDFtransform::KindOf(wcs.Tag()) == "gwcs/wcs" && shiftGWCS(wcs, origin, size))
                    shifted.push_back(scope.empty() ? path : scope + "/" + path);
                else
                    doomed.push_back(path);
            }

            // Last first, so that the indices of earlier sequence items still hold
            for(auto it = doomed.rbegin(); it != doomed.rend(); ++it)
                top->ErasePath(*it);

            std::vector<std::string> removed;
            for(const std::string& path : doomed)
                removed.push_back(scope.empty() ? path : scope + "/" + path);
            return removed;
        }

        /**
         * @brief Remove the GWCS objects found below a node.
         * @param tree  Root of the tree.
         * @param scope Path of the node to search; empty for the whole tree.
         * @return Paths of the objects removed.
         */
        std::vector<std::string> removeGWCS(ASDFnode& tree, const std::string& scope)
        {
            std::vector<std::string> none;
            return fitGWCS(tree, scope, {}, {}, none);
        }

        /**
         * @brief Mask the pixels flagged in a data quality image, if it holds type @p T.
         * @tparam T Pixel type to try; integers only.
         * @param dq    The data quality image.
         * @param bits  Flags that mask a pixel.
         * @param mask  Receives true for each flagged pixel.
         * @param count Receives the number of flagged pixels.
         * @return False if the image holds another type.
         */
        template<typename T>
        bool flagsOf(const FITScube& dq, uint64_t bits, std::valarray<bool>& mask, size_t& count)
        {
            const std::valarray<T>* flags = dq.GetData<T>();
            if(!flags)
                return false;
            mask.resize(flags->size(), false);
            count = 0;
            for(size_t i = 0; i < flags->size(); i++)
                if(static_cast<uint64_t>((*flags)[i]) & bits)
                {
                    mask[i] = true;
                    count++;
                }
            return true;
        }
    }

    namespace asdfdetail
    {
        /** @brief The GWCS of an image, read when first asked for. */
        struct GwcsCache
        {
            std::mutex lock;                        //!< Guards the rest: the GWCS is asked for from several threads
            bool known = false;                     //!< The tree has been searched
            std::string path;                       //!< Path of the GWCS object; empty if none
            std::shared_ptr<const ASDFgwcs> gwcs;   //!< The GWCS; null if none or unusable

            /** @brief Take what another cache knows. @param other The other cache; its lock is held. */
            void copyFrom(const GwcsCache& other)
            {
                known = other.known;
                path = other.path;
                gwcs = other.gwcs;  // Immutable: sharing it is safe
            }
        };
    }

#pragma region - ASDFcube
#pragma region * ctor
    ASDFcube::ASDFcube(): flink(std::make_unique<HeaderLink>()), fgwcs(std::make_unique<GwcsCache>()) {}

    ASDFcube::~ASDFcube() = default;

    ASDFcube::ASDFcube(const ASDFcube& other):
        ftree(other.ftree.Clone()),
        fimagePath(other.fimagePath),
        fcompanions(other.fcompanions),
        flink(std::make_unique<HeaderLink>(*other.flink)),
        fsnapshot(other.fsnapshot),
        fhistory(other.fhistory),
        faligned(other.faligned),
        fgwcs(std::make_unique<GwcsCache>())
    {
        std::lock_guard<std::mutex> guard(other.fgwcs->lock);
        fgwcs->copyFrom(*other.fgwcs);
    }

    ASDFcube& ASDFcube::operator=(const ASDFcube& other)
    {
        if(this == &other)
            return *this;
        ftree = other.ftree.Clone();
        fimagePath = other.fimagePath;
        fcompanions = other.fcompanions;
        *flink = *other.flink;
        FITShdu snapshot(other.fsnapshot);
        fsnapshot.swap(snapshot);
        fhistory = other.fhistory;
        faligned = other.faligned;
        std::scoped_lock guard(fgwcs->lock, other.fgwcs->lock);
        fgwcs->copyFrom(*other.fgwcs);
        return *this;
    }

#pragma endregion
#pragma region * GWCS
    void ASDFcube::forgetGWCS() const
    {
        std::lock_guard<std::mutex> guard(fgwcs->lock);
        fgwcs->known = false;
        fgwcs->path.clear();
        fgwcs->gwcs.reset();
    }

    std::shared_ptr<const ASDFgwcs> ASDFcube::GWCS() const
    {
        std::lock_guard<std::mutex> guard(fgwcs->lock);
        if(fgwcs->known)
            return fgwcs->gwcs;
        fgwcs->known = true;

        const std::string scope = fimagePath.empty() ? std::string() : ASDFnode::ParentPath(fimagePath);
        const ASDFnode* top = ftree.FindPath(scope);
        if(!top)
            return nullptr;
        for(const std::string& path : gwcsObjects(*top))
            if(ASDFtransform::KindOf(top->At(path).Tag()) == "gwcs/wcs")
            {
                fgwcs->path = scope.empty() ? path : scope + "/" + path;
                break;
            }
        if(fgwcs->path.empty())
            return nullptr;

        try
        {
            std::shared_ptr<const ASDFgwcs> gwcs = ASDFgwcs::Read(ftree.At(fgwcs->path));
            if(gwcs->PixelAxes() != Cube().GetDimension())
                std::cerr << "\033[33m[WARNING]\033[0m ASDFcube::GWCS: " << fgwcs->path << " has " << gwcs->PixelAxes()
                          << " pixel axes, the image " << Cube().GetDimension() << "; the FITS WCS of the header is used instead" << std::endl;
            else
                fgwcs->gwcs = std::move(gwcs);
        }
        catch(const FITSexception& e)
        {
            std::cerr << "\033[33m[WARNING]\033[0m ASDFcube::GWCS: " << fgwcs->path << " can't be used: " << e.what()
                      << "; the FITS WCS of the header is used instead" << std::endl;
        }
        return fgwcs->gwcs;
    }

    std::string ASDFcube::GWCSPath() const
    {
        GWCS();
        std::lock_guard<std::mutex> guard(fgwcs->lock);
        return fgwcs->path;
    }

#pragma endregion
#pragma region * Tree
    void ASDFcube::Set(const std::string& path, ASDFnode value)
    {
        SyncFromHDU();
        ftree.SetPath(path, std::move(value));
        RefreshHDU();
    }

    bool ASDFcube::Erase(const std::string& path)
    {
        SyncFromHDU();
        const bool erased = ftree.ErasePath(path);
        RefreshHDU();
        return erased;
    }

#pragma endregion
#pragma region * Companions
    std::vector<size_t> ASDFcube::axesNow() const
    {
        return axesOf(Cube());
    }

    std::vector<std::pair<std::string, std::shared_ptr<ASDFarray>>> ASDFcube::companions(const ASDFnode& tree) const
    {
        std::vector<std::pair<std::string, std::shared_ptr<ASDFarray>>> found;
        if(!fcompanions)
            return found;

        const std::string parent = ASDFnode::ParentPath(fimagePath);
        const ASDFnode* mapping = tree.FindPath(parent);
        if(!mapping || !mapping->IsMapping())
            return found;

        for(size_t i = 0; i < mapping->Size(); i++)
        {
            const ASDFnode& child = mapping->Items()[i];
            const std::string path = ASDFnode::ChildPath(parent, mapping->Keys()[i]);
            if(child.IsArray() && path != fimagePath)
                found.emplace_back(path, child.Array());
        }
        return found;
    }

    std::vector<std::string> ASDFcube::CompanionNames() const
    {
        // Companions left behind by a change of shape the hooks didn't see are no longer the image's
        const bool stale = faligned != axesNow();
        std::vector<std::string> names;
        for(const auto& [path, array] : companions(ftree))
            if(!stale || array->Axes() != faligned)
                names.push_back(ASDFnode::SplitPath(path).back());
        return names;
    }

    std::shared_ptr<ASDFarray> ASDFcube::Companion(const std::string& name) const
    {
        for(const std::string& known : CompanionNames())
            if(known == name)
                return ftree.At(ASDFnode::ChildPath(ASDFnode::ParentPath(fimagePath), name)).Array();
        return nullptr;
    }

    size_t ASDFcube::MaskFromDQ(uint64_t bits, const std::string& name)
    {
        std::shared_ptr<ASDFarray> dq = Companion(name);
        if(!dq)
            throw FITSexception(KEY_NO_EXIST, kWho, "MaskFromDQ", "the image has no companion '"+name+"'");
        if(dq->Axes() != axesNow())
            throw FITSexception(BAD_DIMEN, kWho, "MaskFromDQ", "'"+name+"' is "+axesText(dq->Axes())+", the image "+axesText(axesNow()));

        const bool wasLoaded = dq->IsLoaded();
        std::shared_ptr<FITScube> flags = dq->Load();

        std::valarray<bool> mask;
        size_t count = 0;
        const bool done = flagsOf<uint8_t>(*flags, bits, mask, count)  || flagsOf<uint16_t>(*flags, bits, mask, count)
                       || flagsOf<uint32_t>(*flags, bits, mask, count) || flagsOf<uint64_t>(*flags, bits, mask, count)
                       || flagsOf<int8_t>(*flags, bits, mask, count)   || flagsOf<int16_t>(*flags, bits, mask, count)
                       || flagsOf<int32_t>(*flags, bits, mask, count)  || flagsOf<int64_t>(*flags, bits, mask, count);
        flags.reset();
        if(!wasLoaded)
            dq->Release();
        if(!done)
            throw FITSexception(BAD_DATATYPE, kWho, "MaskFromDQ", "'"+name+"' doesn't hold integers");

        if(count > 0)
            Cube().MaskPixels(mask);
        return count;
    }

#pragma endregion
#pragma region * Header and tree
    void ASDFcube::SyncFromHDU()
    {
        syncHeader(ftree, *flink, Cube().HDU(), fsnapshot);
        FITShdu snapshot(Cube().HDU());
        fsnapshot.swap(snapshot);
    }

    void ASDFcube::RefreshHDU()
    {
        forgetGWCS();   // Meta() may have been changed
        SyncFromHDU();
        FITShdu& hdu = Cube().HDU();
        rebuildHeader(ftree, *flink, hdu);
        FITShdu snapshot(hdu);
        fsnapshot.swap(snapshot);
        if(hdu.Exists("CTYPE1") || hdu.Exists("CRPIX1") || hdu.Exists("WCSAXES"))
            Cube().reLoadWCS();
    }

#pragma endregion
#pragma region * Changes of shape
    void ASDFcube::alignTo(const std::vector<size_t>& before)
    {
        if(faligned == before)
            return;

        // The shape changed without the hooks: what matched the old shape no longer matches
        forgetGWCS();
        std::vector<std::string> dropped;
        for(const auto& [path, array] : companions(ftree))
            if(array->Axes() == faligned)
            {
                ftree.ErasePath(path);
                dropped.push_back(ASDFnode::SplitPath(path).back());
            }
        const std::vector<std::string> gwcs = removeGWCS(ftree, ASDFnode::ParentPath(fimagePath));

        std::string note = "Image reshaped from "+axesText(faligned)+" to "+axesText(before);
        if(!dropped.empty())
            note += "; companions dropped: "+join(dropped);
        if(!gwcs.empty())
            note += "; GWCS removed: "+join(gwcs);
        fhistory.push_back(note);
        faligned = before;
    }

    void ASDFcube::AfterWindow(const std::vector<size_t>& before, const std::vector<size_t>& origin, const std::string& what)
    {
        const std::vector<size_t> now = axesNow();
        if(now == before)
            return;
        alignTo(before);

        std::vector<size_t> start(origin);
        start.resize(before.size(), 0);

        std::vector<std::string> cut;
        for(const auto& [path, array] : companions(ftree))
            if(array->Axes() == before)
            {
                ftree.At(path) = ASDFnode(array->Window(start, now));
                cut.push_back(ASDFnode::SplitPath(path).back());
            }

        // A 2D GWCS follows the window; any other is removed
        forgetGWCS();
        std::vector<std::string> shifted;
        const std::vector<std::string> gwcs = before.size() == 2
            ? fitGWCS(ftree, ASDFnode::ParentPath(fimagePath), start, now, shifted)
            : removeGWCS(ftree, ASDFnode::ParentPath(fimagePath));

        std::string note = what+": "+axesText(before)+" to "+axesText(now)+" from pixel "+axesText(start);
        if(!cut.empty())
            note += "; companions cut too: "+join(cut);
        if(!shifted.empty())
            note += "; GWCS shifted to the window: "+join(shifted);
        if(!gwcs.empty())
            note += "; GWCS removed: "+join(gwcs);
        fhistory.push_back(note);
        faligned = now;
    }

    void ASDFcube::AfterLayer(const std::vector<size_t>& before, size_t index, const std::string& what)
    {
        const std::vector<size_t> now = axesNow();
        alignTo(before);

        std::vector<std::string> cut;
        for(const auto& [path, array] : companions(ftree))
            if(array->Axes() == before)
            {
                ftree.At(path) = ASDFnode(array->Layer(index));
                cut.push_back(ASDFnode::SplitPath(path).back());
            }
        forgetGWCS();
        const std::vector<std::string> gwcs = removeGWCS(ftree, ASDFnode::ParentPath(fimagePath));

        std::string note = what+": layer "+std::to_string(index)+" of "+axesText(before);
        if(!cut.empty())
            note += "; companions cut too: "+join(cut);
        if(!gwcs.empty())
            note += "; GWCS removed: "+join(gwcs);
        fhistory.push_back(note);
        faligned = now;
    }

    void ASDFcube::AfterReshape(const std::vector<size_t>& before, const std::string& what)
    {
        const std::vector<size_t> now = axesNow();
        if(now == before)
            return;
        alignTo(before);

        std::vector<std::string> dropped;
        for(const auto& [path, array] : companions(ftree))
            if(array->Axes() == before)
            {
                ftree.ErasePath(path);
                dropped.push_back(ASDFnode::SplitPath(path).back());
            }
        forgetGWCS();
        const std::vector<std::string> gwcs = removeGWCS(ftree, ASDFnode::ParentPath(fimagePath));

        std::string note = what+": "+axesText(before)+" to "+axesText(now);
        if(!dropped.empty())
            note += "; companions dropped: "+join(dropped);
        if(!gwcs.empty())
            note += "; GWCS removed: "+join(gwcs);
        fhistory.push_back(note);
        faligned = now;
    }

    void ASDFcube::AfterProcessing(const std::vector<size_t>& before, const std::string& what)
    {
        const std::vector<size_t> now = axesNow();
        alignTo(before);

        std::vector<std::string> dropped;
        for(const auto& [path, array] : companions(ftree))
            if(array->Axes() == before)
            {
                ftree.ErasePath(path);
                dropped.push_back(ASDFnode::SplitPath(path).back());
            }
        forgetGWCS();
        std::vector<std::string> gwcs;
        if(now != before)
            gwcs = removeGWCS(ftree, ASDFnode::ParentPath(fimagePath));

        std::string note = what;
        if(now != before)
            note += ": "+axesText(before)+" to "+axesText(now);
        if(!dropped.empty())
            note += "; companions dropped: "+join(dropped);
        if(!gwcs.empty())
            note += "; GWCS removed: "+join(gwcs);
        fhistory.push_back(note);
        faligned = now;
    }

#pragma endregion
#pragma region * History and writing
    void ASDFcube::AddHistory(const std::string& text)
    {
        fhistory.push_back(text);
    }

    void ASDFcube::WriteASDF(const std::string& path, const ASDFwriteOptions& options) const
    {
        checkCompression(options.compression);

        // Work on copies: writing changes nothing of this object
        ASDFnode tree = ftree;
        HeaderLink link = *flink;
        std::vector<std::string> notes = fhistory;

        const std::vector<size_t> now = axesNow();
        if(faligned != now)
        {
            std::vector<std::string> dropped;
            for(const auto& [companionPath, array] : companions(tree))
                if(array->Axes() == faligned)
                {
                    tree.ErasePath(companionPath);
                    dropped.push_back(ASDFnode::SplitPath(companionPath).back());
                }
            const std::vector<std::string> gwcs = removeGWCS(tree, ASDFnode::ParentPath(fimagePath));
            std::string note = "Image reshaped from "+axesText(faligned)+" to "+axesText(now);
            if(!dropped.empty())
                note += "; companions dropped: "+join(dropped);
            if(!gwcs.empty())
                note += "; GWCS removed: "+join(gwcs);
            notes.push_back(note);
        }

        PixelBlock pixels = cubePixels(Cube());
        FITShdu header(Cube().HDU());
        if(!link.cardsPath.empty())
        {
            // The astropy layout keeps the whole FITS header: its structure follows the pixels written
            setImageStructure(header, pixels.type, std::vector<uint64_t>(now.begin(), now.end()));
            if(header.Exists("BLANK"))
                header.DeleteKey("BLANK");
            if(pixels.blank.GetKind() == ASDFnode::Kind::unsigned_integer)
                header.ValueForKey("BLANK", pixels.blank.Get<uint64_t>());
            else if(pixels.blank.IsInteger())
                header.ValueForKey("BLANK", pixels.blank.Get<int64_t>());

            syncHeader(tree, link, header, fsnapshot);
            if(ASDFnode* cards = tree.FindPath(link.cardsPath))
                for(const char* stale : {"CHECKSUM", "DATASUM"})
                    removeCards(*cards, stale);
        }
        else
        {
            syncHeader(tree, link, header, fsnapshot);

            // Masked pixels of an integer image hold BLANK: say which value, so that they are masked again on reading
            auto known = link.keyPaths.find("BLANK");
            if(!pixels.blank.IsNull())
                tree.SetPath(known != link.keyPaths.end() ? known->second : ASDFnode::ChildPath(link.sinkPath, "BLANK"), pixels.blank);
            else if(known != link.keyPaths.end())
                tree.ErasePath(known->second);
        }

        for(const std::string& note : notes)
            appendHistory(tree, historyEntry(note));
        if(!options.history.empty())
            appendHistory(tree, historyEntry(options.history));

        writeTreeFile(tree, fimagePath, std::move(pixels), options.compression, path, options.replace);
    }
#pragma endregion
#pragma endregion

#pragma region - ASDFimg<T>
    template<typename T>
    ASDFimg<T>::ASDFimg(const std::vector<size_t>& axes): FITSimg<T>(axes), ASDFcube() {}

    template<typename T>
    ASDFimg<T>::ASDFimg(const ASDFimg<T>& other): FITSimg<T>(other), ASDFcube(other) {}

    template<typename T>
    ASDFimg<T>& ASDFimg<T>::operator=(const ASDFimg<T>& other)
    {
        if(this != &other)
        {
            FITSimg<T>::operator=(other);
            ASDFcube::operator=(other);
        }
        return *this;
    }

    template<typename T>
    ASDFimg<T>::~ASDFimg() = default;

    template<typename T>
    std::shared_ptr<ASDFimg<T>> ASDFimg<T>::Adopt(FITSimg<T>& image, const ASDFcube& model, bool inheritKeywords)
    {
        // A one-pixel image receives the pixels: no copy of them is made
        auto out = std::make_shared<ASDFimg<T>>(std::vector<size_t>{1});
        std::valarray<T>* mine = out->template GetData<T>();
        std::valarray<T>* theirs = image.template GetData<T>();
        if(mine && theirs)
            std::swap(*mine, *theirs);

        out->Naxis = axesOf(image);
        out->mask = image.GetMask();
        out->HDU().swap(image.HDU());

        // A header built afresh lacks the keywords read from the tree: they weren't removed on purpose
        if(inheritKeywords)
        {
            const FITShdu& from = model.Cube().HDU();
            for(auto it = from.begin(); it != from.end(); ++it)
                if(!out->HDU().Exists(it->first) && !isStructuralKey(it->first, false))
                    out->HDU().ValueForKey(it->first, it->second.value(), it->second.type(), it->second.comment());
        }
        if(out->HDU().Exists("EXTNAME"))
            out->name = out->HDU().GetValueForKey("EXTNAME");
        if constexpr (!std::is_floating_point_v<T>)
        {
            if(out->HDU().Exists("BLANK"))
                out->Blank(std::is_unsigned_v<T> ? static_cast<T>(out->HDU().GetUInt64ValueForKey("BLANK"))
                                                 : static_cast<T>(out->HDU().GetInt64ValueForKey("BLANK")));
        }
        out->reLoadWCS();
        out->ASDFcube::operator=(model);
        return out;
    }

    template<typename T>
    std::shared_ptr<FITScube> ASDFimg<T>::Layer(const size_t& index) const
    {
        std::shared_ptr<FITScube> plain = FITSimg<T>::Layer(index);
        std::shared_ptr<ASDFimg<T>> out = Adopt(dynamic_cast<FITSimg<T>&>(*plain), *this);
        out->AfterLayer(axesOf(*this), index, "Layer");
        return out;
    }

    template<typename T>
    std::shared_ptr<FITScube> ASDFimg<T>::Window(size_t xMin, size_t yMin, size_t xSize, size_t ySize) const
    {
        std::shared_ptr<FITScube> plain = FITSimg<T>::Window(xMin, yMin, xSize, ySize);
        std::shared_ptr<ASDFimg<T>> out = Adopt(dynamic_cast<FITSimg<T>&>(*plain), *this);
        out->AfterWindow(axesOf(*this), {xMin, yMin}, "Window");
        return out;
    }

    template<typename T>
    std::shared_ptr<FITScube> ASDFimg<T>::Rebin(const std::vector<size_t>& nbin, bool doMean) const
    {
        // FITSimg gives the result a header of its own, with the structure and the WCS only
        std::shared_ptr<FITScube> plain = FITSimg<T>::Rebin(nbin, doMean);
        std::shared_ptr<ASDFimg<T>> out = Adopt(dynamic_cast<FITSimg<T>&>(*plain), *this, true);
        out->AfterReshape(axesOf(*this), doMean ? "Rebin (mean)" : "Rebin (sum)");
        return out;
    }

    template<typename T>
    std::shared_ptr<FITScube> ASDFimg<T>::Overlay(const FITScube::overlay& method, const std::pair<double,double>& clip) const
    {
        std::shared_ptr<FITScube> plain = FITSimg<T>::Overlay(method, clip);
        std::shared_ptr<ASDFimg<T>> out = Adopt(dynamic_cast<FITSimg<T>&>(*plain), *this);
        out->AfterReshape(axesOf(*this), "Overlay of the layers");
        return out;
    }

    template<typename T>
    void ASDFimg<T>::Resize(const size_t& xMin, const size_t& yMin, const size_t& xSize, const size_t& ySize)
    {
        const std::vector<size_t> before = axesOf(*this);
        FITSimg<T>::Resize(xMin, yMin, xSize, ySize);
        AfterWindow(before, {xMin, yMin}, "Resize");
    }

    template<typename T>
    void ASDFimg<T>::AddLayer(const FITScube& layer)
    {
        const std::vector<size_t> before = axesOf(*this);
        FITSimg<T>::AddLayer(layer);
        AfterReshape(before, "AddLayer");
    }

    template<typename T>
    void ASDFimg<T>::AddLayer(const FITSimg<T>& layer)
    {
        const std::vector<size_t> before = axesOf(*this);
        FITSimg<T>::AddLayer(layer);
        AfterReshape(before, "AddLayer");
    }

    template<typename T>
    worldCoords ASDFimg<T>::WorldCoordinates(const pixelCoords& pixel, const int& wcsIndex) const
    {
        if(wcsIndex == 0)
            if(std::shared_ptr<const ASDFgwcs> gwcs = GWCS())
                return gwcs->PixelToWorld(pixel);
        return FITSimg<T>::WorldCoordinates(pixel, wcsIndex);
    }

    template<typename T>
    worldVectors ASDFimg<T>::WorldCoordinatesVector(const pixelVectors& pixels, const int& wcsIndex) const
    {
        if(wcsIndex == 0)
            if(std::shared_ptr<const ASDFgwcs> gwcs = GWCS())
            {
                worldVectors world;
                world.reserve(pixels.size());
                for(const pixelCoords& pixel : pixels)
                    world.push_back(gwcs->PixelToWorld(pixel));
                return world;
            }
        return FITSimg<T>::WorldCoordinatesVector(pixels, wcsIndex);
    }

    template<typename T>
    pixelCoords ASDFimg<T>::World2Pixel(const worldCoords& world, const int& wcsIndex) const
    {
        if(wcsIndex == 0)
            if(std::shared_ptr<const ASDFgwcs> gwcs = GWCS())
                return gwcs->WorldToPixel(world);
        return FITSimg<T>::World2Pixel(world, wcsIndex);
    }

    template<typename T>
    pixelVectors ASDFimg<T>::World2PixelVector(const worldVectors& world, const int& wcsIndex) const
    {
        if(wcsIndex == 0)
            if(std::shared_ptr<const ASDFgwcs> gwcs = GWCS())
            {
                pixelVectors pixels;
                pixels.reserve(world.size());
                for(const worldCoords& position : world)
                    pixels.push_back(gwcs->WorldToPixel(position));
                return pixels;
            }
        return FITSimg<T>::World2PixelVector(world, wcsIndex);
    }

    template<typename T>
    bool ASDFimg<T>::HasWorldCoordinates() const
    {
        return GWCS() != nullptr || FITSimg<T>::HasWorldCoordinates();
    }

    template<typename T>
    std::vector<std::string> ASDFimg<T>::WorldAxisNames(const int& wcsIndex) const
    {
        if(wcsIndex == 0)
            if(std::shared_ptr<const ASDFgwcs> gwcs = GWCS())
                return gwcs->OutputFrame().axesNames;
        return FITSimg<T>::WorldAxisNames(wcsIndex);
    }

    template<typename T>
    std::vector<std::string> ASDFimg<T>::WorldAxisUnits(const int& wcsIndex) const
    {
        if(wcsIndex == 0)
            if(std::shared_ptr<const ASDFgwcs> gwcs = GWCS())
                return gwcs->OutputFrame().units;
        return FITSimg<T>::WorldAxisUnits(wcsIndex);
    }

    template class ASDFimg<int8_t>;
    template class ASDFimg<uint8_t>;
    template class ASDFimg<int16_t>;
    template class ASDFimg<uint16_t>;
    template class ASDFimg<int32_t>;
    template class ASDFimg<uint32_t>;
    template class ASDFimg<int64_t>;
    template class ASDFimg<uint64_t>;
    template class ASDFimg<float>;
    template class ASDFimg<double>;

#pragma endregion
#pragma region - asdfLike and the access for the library
    namespace
    {
        /**
         * @brief Give a result the tree of a model, if the result holds pixels of type @p T.
         * @tparam T Pixel type to try.
         * @param result The result.
         * @param model  The model.
         * @param what   Description of how the result was made.
         * @param out    Receives the ASDF image.
         * @return False if the result holds another type.
         */
        template<typename T>
        bool likeAs(std::shared_ptr<FITScube>& result, const ASDFcube& model, const std::string& what, std::shared_ptr<FITScube>& out)
        {
            FITSimg<T>* image = dynamic_cast<FITSimg<T>*>(result.get());
            if(!image)
                return false;

            std::shared_ptr<ASDFimg<T>> adopted;
            if(result.use_count() == 1)
                adopted = ASDFimg<T>::Adopt(*image, model, true);   // nothing else sees it: take its pixels
            else
            {
                FITSimg<T> copy(*image);
                adopted = ASDFimg<T>::Adopt(copy, model, true);
            }
            CubeAccess::afterProcessing(*adopted, axesOf(model.Cube()), what);
            out = adopted;
            return true;
        }
    }

    std::shared_ptr<FITScube> asdfLike(std::shared_ptr<FITScube> result, const ASDFcube& model, const std::string& what)
    {
        if(!result)
            throw FITSexception(SHARED_NULPTR, kWho, "asdfLike", "no result given");
        if(dynamic_cast<ASDFcube*>(result.get()))
            return result;

        std::shared_ptr<FITScube> out;
        const bool done = likeAs<int8_t>(result, model, what, out)   || likeAs<uint8_t>(result, model, what, out)
                       || likeAs<int16_t>(result, model, what, out)  || likeAs<uint16_t>(result, model, what, out)
                       || likeAs<int32_t>(result, model, what, out)  || likeAs<uint32_t>(result, model, what, out)
                       || likeAs<int64_t>(result, model, what, out)  || likeAs<uint64_t>(result, model, what, out)
                       || likeAs<float>(result, model, what, out)    || likeAs<double>(result, model, what, out);
        if(!done)
            throw FITSexception(BAD_DATATYPE, kWho, "asdfLike", "the result holds no pixels of a known type");
        return out;
    }

    void WriteASDF(const ASDFnode& tree, const std::string& path, const ASDFwriteOptions& options)
    {
        checkCompression(options.compression);
        ASDFnode copy = tree;
        if(!options.history.empty())
            appendHistory(copy, historyEntry(options.history));
        writeTreeFile(copy, std::string(), PixelBlock(), options.compression, path, options.replace);
    }

    namespace asdfdetail
    {
        void CubeAccess::attach(ASDFcube& cube, ImageState state)
        {
            cube.ftree = std::move(state.tree);
            cube.fimagePath = std::move(state.imagePath);
            cube.fcompanions = state.companions;
            *cube.flink = std::move(state.link);

            HeaderLink& link = *cube.flink;
            if(link.cardsPath.empty())
            {
                FITShdu scratch;
                KeyPaths paths;
                if(link.hasKeywords)
                    if(const ASDFnode* mapping = cube.ftree.FindPath(link.keywordsPath))
                        readTreeKeywords(*mapping, std::string(), link.keywordsPath, scratch, &paths);
                if(link.rootKeywords)
                    readRootKeywords(cube.ftree, scratch, &paths);
                link.keyPaths = std::move(paths);
            }

            FITShdu snapshot(cube.Cube().HDU());
            cube.fsnapshot.swap(snapshot);
            cube.faligned = cube.axesNow();
            cube.forgetGWCS();
        }
    }
#pragma endregion
}
