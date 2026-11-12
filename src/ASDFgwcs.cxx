//
//  ASDFgwcs.cxx
//  asdf2fits
//
//  Written by Claude (Anthropic's AI model) for William Gillard.
//  Copyright (c) 2026 William Gillard.
//
//  SPDX-License-Identifier: BSD-3-Clause
//  Licence text: see the LICENSE file at the root of this project.

/**
 * @file ASDFgwcs.cxx
 * @brief Implementation of DSL::ASDFgwcs.
 */
#include "ASDFgwcs.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <limits>

#include <DSTfits/FITSimg.h>

namespace DSL
{
    namespace
    {
        const char* const kWho = "ASDFgwcs";    //!< Name given in the exceptions thrown here
        const double kNaN = std::numeric_limits<double>::quiet_NaN();

        /**
         * @brief Read a list of texts.
         * @param node The list; may be null.
         * @return The texts; empty if @p node isn't a list.
         */
        std::vector<std::string> texts(const ASDFnode* node)
        {
            std::vector<std::string> out;
            if(!node || !node->IsSequence())
                return out;
            for(const ASDFnode& item : node->Items())
                out.push_back(item.IsString() ? item.Get<std::string>() : item.Text());
            return out;
        }

        /**
         * @brief Read a frame.
         * @param node Node of the frame: a mapping, or only a name in old files.
         * @return The frame.
         */
        ASDFframe readFrame(const ASDFnode& node)
        {
            ASDFframe frame;
            if(node.IsString())
            {
                frame.name = node.Get<std::string>();
                return frame;
            }
            if(!node.IsMapping())
                return frame;

            frame.kind = ASDFtransform::KindOf(node.Tag());
            if(const ASDFnode* name = node.Find("name"); name && name->IsString())
                frame.name = name->Get<std::string>();

            if(const ASDFnode* parts = node.Find("frames"); parts && parts->IsSequence())
            {
                // A composite frame: the axes of its frames, one after the other
                for(const ASDFnode& part : parts->Items())
                {
                    const ASDFframe sub = readFrame(part);
                    frame.axesNames.insert(frame.axesNames.end(), sub.axesNames.begin(), sub.axesNames.end());
                    frame.units.insert(frame.units.end(), sub.units.begin(), sub.units.end());
                    if(frame.referenceFrame.empty())
                        frame.referenceFrame = sub.referenceFrame;
                }
            }
            else
            {
                frame.axesNames = texts(node.Find("axes_names"));
                frame.units = texts(node.Find("unit"));
            }

            if(const ASDFnode* reference = node.Find("reference_frame"); reference && !reference->Tag().empty())
            {
                const std::string kind = ASDFtransform::KindOf(reference->Tag());
                frame.referenceFrame = kind.substr(kind.rfind('/') + 1);
            }
            return frame;
        }

        /**
         * @brief Read the bounding box of the first transform.
         * @details Two forms: a bounding_box object whose intervals are named after the inputs of the transform,
         * or, in older files, a list of [lower, upper] in the C order, the last input first.
         * @param transform Node of the first transform.
         * @param inputs    Its number of inputs.
         * @return The range of each input, first input first; empty if there is no bounding box.
         */
        std::vector<std::pair<double,double>> readBox(const ASDFnode& transform, size_t inputs)
        {
            const ASDFnode* box = transform.Find("bounding_box");
            if(!box || box->IsNull())
                return {};

            const double inf = std::numeric_limits<double>::infinity();
            std::vector<std::pair<double,double>> out(inputs, {-inf, inf});

            const auto range = [](const ASDFnode& node)
            {
                const std::vector<double> values = ASDFtransform::Numbers(node);
                if(values.size() != 2)
                    throw WCSexception(BAD_OPTION, kWho, "readBox", "an interval needs two values");
                return std::make_pair(std::min(values[0], values[1]), std::max(values[0], values[1]));
            };

            if(box->IsMapping())
            {
                const ASDFnode* intervals = box->Find("intervals");
                if(!intervals || !intervals->IsMapping())
                    throw WCSexception(BAD_OPTION, kWho, "readBox", "the bounding box has no intervals");
                const std::vector<std::string> names = texts(transform.Find("inputs"));
                for(size_t k = 0; k < intervals->Size(); k++)
                {
                    const std::string& key = intervals->Keys()[k];
                    size_t index = k;
                    auto it = std::find(names.begin(), names.end(), key);
                    if(it != names.end())
                        index = static_cast<size_t>(it - names.begin());
                    else if(key.size() > 1 && key[0] == 'x' && std::all_of(key.begin() + 1, key.end(), [](char c) { return std::isdigit(static_cast<unsigned char>(c)); }))
                        index = std::stoul(key.substr(1));
                    if(index >= inputs)
                        throw WCSexception(BAD_OPTION, kWho, "readBox", "the bounding box names an unknown input '"+key+"'");
                    out[index] = range(intervals->Items()[k]);
                }
                return out;
            }

            const std::vector<double> values = ASDFtransform::Numbers(*box);
            if(values.size() != 2*inputs)
                throw WCSexception(BAD_OPTION, kWho, "readBox", "the bounding box needs two values per input");
            for(size_t k = 0; k < inputs; k++)
            {
                const size_t index = inputs - 1 - k;
                out[index] = std::make_pair(std::min(values[2*k], values[2*k + 1]), std::max(values[2*k], values[2*k + 1]));
            }
            return out;
        }

        /** @brief Bring a difference of longitudes into [-180, 180). @param d The difference, degrees. @return The same angle. */
        double wrap180(double d)
        {
            d = std::fmod(d + 180., 360.);
            if(d < 0.)
                d += 360.;
            return d - 180.;
        }

        /**
         * @brief Solve a linear system by Gaussian elimination with partial pivoting.
         * @param a Matrix, n x n, C order; destroyed.
         * @param b Right-hand side; receives the solution.
         * @return False if the matrix is singular.
         */
        bool solveLinear(std::vector<double>& a, std::vector<double>& b)
        {
            const size_t n = b.size();
            for(size_t col = 0; col < n; col++)
            {
                size_t pivot = col;
                for(size_t r = col + 1; r < n; r++)
                    if(std::fabs(a[r*n + col]) > std::fabs(a[pivot*n + col]))
                        pivot = r;
                if(!(std::fabs(a[pivot*n + col]) > 0.))
                    return false;
                for(size_t c = 0; c < n; c++)
                    std::swap(a[col*n + c], a[pivot*n + c]);
                std::swap(b[col], b[pivot]);
                for(size_t r = col + 1; r < n; r++)
                {
                    const double factor = a[r*n + col] / a[col*n + col];
                    for(size_t c = col; c < n; c++)
                        a[r*n + c] -= factor * a[col*n + c];
                    b[r] -= factor * b[col];
                }
            }
            for(size_t r = n; r-- > 0;)
            {
                double sum = b[r];
                for(size_t c = r + 1; c < n; c++)
                    sum -= a[r*n + c] * b[c];
                b[r] = sum / a[r*n + r];
            }
            return true;
        }
    }

    std::shared_ptr<const ASDFgwcs> ASDFgwcs::Read(const ASDFnode& wcs)
    {
        if(!wcs.IsMapping())
            throw WCSexception(BAD_OPTION, kWho, "Read", "a GWCS must be a mapping");
        const ASDFnode* steps = wcs.Find("steps");
        if(!steps || !steps->IsSequence() || steps->Size() < 2)
            throw WCSexception(BAD_OPTION, kWho, "Read", "a GWCS needs 'steps', with two frames at least");

        std::shared_ptr<ASDFgwcs> gwcs(new ASDFgwcs());
        if(const ASDFnode* name = wcs.Find("name"); name && name->IsString())
            gwcs->fname = name->Get<std::string>();

        const size_t last = steps->Size() - 1;
        for(size_t k = 0; k <= last; k++)
        {
            const ASDFnode& step = steps->Items()[k];
            if(!step.IsMapping())
                throw WCSexception(BAD_OPTION, kWho, "Read", "step "+std::to_string(k)+" isn't a mapping");
            const ASDFnode* frame = step.Find("frame");
            gwcs->fframes.push_back(frame ? readFrame(*frame) : ASDFframe());
            if(k == last)
                break;  // The transform of the last step is null: the world goes nowhere

            const ASDFnode* transform = step.Find("transform");
            if(!transform || transform->IsNull())
                throw WCSexception(BAD_OPTION, kWho, "Read", "step "+std::to_string(k)+" ("+gwcs->fframes.back().name+") has no transform");
            try
            {
                gwcs->fsteps.push_back(ASDFtransform::Read(*transform));
            }
            catch(const FITSexception& e)
            {
                throw WCSexception(BAD_OPTION, kWho, "Read", "step "+std::to_string(k)+" ("+gwcs->fframes.back().name+"): "+e.what());
            }
            if(k > 0 && gwcs->fsteps[k-1]->Outputs() != gwcs->fsteps[k]->Inputs())
                throw WCSexception(BAD_OPTION, kWho, "Read", "step "+std::to_string(k-1)+" gives "+std::to_string(gwcs->fsteps[k-1]->Outputs())
                                   +" values, step "+std::to_string(k)+" takes "+std::to_string(gwcs->fsteps[k]->Inputs()));
        }

        for(const auto& step : gwcs->fsteps)
        {
            std::shared_ptr<const ASDFtransform> inverse = step->Inverse();
            if(inverse && (inverse->Inputs() != step->Outputs() || inverse->Outputs() != step->Inputs()))
                inverse.reset();
            gwcs->finverses.push_back(std::move(inverse));
        }

        gwcs->fbox = readBox(steps->Items()[0].At("transform"), gwcs->fsteps[0]->Inputs());
        return gwcs;
    }

    size_t ASDFgwcs::PixelAxes() const
    {
        return fsteps.front()->Inputs();
    }

    size_t ASDFgwcs::WorldAxes() const
    {
        return fsteps.back()->Outputs();
    }

    const std::shared_ptr<const ASDFtransform>& ASDFgwcs::Step(size_t step) const
    {
        if(step >= fsteps.size())
            throw WCSexception(BAD_DIMEN, kWho, "Step", "the GWCS has "+std::to_string(fsteps.size())+" steps, no step "+std::to_string(step));
        return fsteps[step];
    }

    bool ASDFgwcs::HasInverse() const
    {
        return std::all_of(finverses.begin(), finverses.end(), [](const auto& inverse) { return inverse != nullptr; });
    }

    bool ASDFgwcs::inside(const std::vector<double>& pixel) const
    {
        for(size_t k = 0; k < fbox.size() && k < pixel.size(); k++)
            if(!(pixel[k] >= fbox[k].first && pixel[k] <= fbox[k].second))
                return false;
        return true;
    }

    size_t ASDFgwcs::frameIndex(const std::string& name) const
    {
        for(size_t k = 0; k < fframes.size(); k++)
            if(fframes[k].name == name)
                return k;
        throw WCSexception(BAD_OPTION, kWho, "frameIndex", "the GWCS has no frame '"+name+"'");
    }

    std::vector<double> ASDFgwcs::run(size_t from, size_t to, std::vector<double> coords) const
    {
        if(from <= to)
        {
            for(size_t k = from; k < to; k++)
            {
                std::vector<double> next(fsteps[k]->Outputs());
                fsteps[k]->Evaluate(coords.data(), next.data());
                coords.swap(next);
            }
            return coords;
        }

        for(size_t k = to; k < from; k++)
            if(!finverses[k])
                return solve(to, from, coords);

        for(size_t k = from; k-- > to;)
        {
            std::vector<double> next(finverses[k]->Outputs());
            finverses[k]->Evaluate(coords.data(), next.data());
            coords.swap(next);
        }
        return coords;
    }

    std::vector<double> ASDFgwcs::solve(size_t from, size_t to, const std::vector<double>& target) const
    {
        // Find x in frame 'from' such that the transforms take it to 'target' in frame 'to'
        const size_t n = fsteps[from]->Inputs();
        std::vector<double> failed(n, kNaN);
        if(fsteps[to - 1]->Outputs() != n)
            return failed;
        for(double value : target)
            if(!std::isfinite(value))
                return failed;

        // Start from the middle of the bounding box, or the origin, taken to frame 'from'
        std::vector<double> x(fsteps[0]->Inputs(), 0.);
        for(size_t k = 0; k < fbox.size() && k < x.size(); k++)
            if(std::isfinite(fbox[k].first) && std::isfinite(fbox[k].second))
                x[k] = 0.5*(fbox[k].first + fbox[k].second);
        x = run(0, from, x);

        const bool celestial = fframes[to].IsCelestial();
        const auto residual = [&](const std::vector<double>& at)
        {
            std::vector<double> r = run(from, to, at);
            for(size_t k = 0; k < n; k++)
                r[k] -= target[k];
            if(celestial)
                r[0] = wrap180(r[0]);
            return r;
        };

        for(int iteration = 0; iteration < 50; iteration++)
        {
            const std::vector<double> r = residual(x);
            for(double value : r)
                if(!std::isfinite(value))
                    return failed;

            // Jacobian by central differences
            std::vector<double> jacobian(n*n);
            for(size_t c = 0; c < n; c++)
            {
                const double h = 1e-6*std::max(1., std::fabs(x[c]));
                std::vector<double> plus(x), minus(x);
                plus[c] += h;
                minus[c] -= h;
                const std::vector<double> rp = residual(plus), rm = residual(minus);
                for(size_t row = 0; row < n; row++)
                {
                    double d = rp[row] - rm[row];
                    if(celestial && row == 0)
                        d = wrap180(d);
                    jacobian[row*n + c] = d / (2.*h);
                }
            }

            std::vector<double> step(r);
            for(double& value : step)
                value = -value;
            if(!solveLinear(jacobian, step))
                return failed;

            bool converged = true;
            for(size_t k = 0; k < n; k++)
            {
                x[k] += step[k];
                if(std::fabs(step[k]) > 1e-11*std::max(1., std::fabs(x[k])))
                    converged = false;
            }
            if(converged)
                return x;
        }
        return failed;
    }

    std::vector<double> ASDFgwcs::PixelToWorld(const std::vector<double>& pixel, bool withBoundingBox) const
    {
        if(pixel.size() != PixelAxes())
            throw WCSexception(BAD_DIMEN, kWho, "PixelToWorld", "the GWCS has "+std::to_string(PixelAxes())+" pixel axes, not "+std::to_string(pixel.size()));
        if(withBoundingBox && !inside(pixel))
            return std::vector<double>(WorldAxes(), kNaN);
        return run(0, fsteps.size(), pixel);
    }

    std::vector<double> ASDFgwcs::WorldToPixel(const std::vector<double>& world, bool withBoundingBox) const
    {
        if(world.size() != WorldAxes())
            throw WCSexception(BAD_DIMEN, kWho, "WorldToPixel", "the GWCS has "+std::to_string(WorldAxes())+" world axes, not "+std::to_string(world.size()));
        std::vector<double> pixel = run(fsteps.size(), 0, world);
        if(withBoundingBox && !inside(pixel))
            std::fill(pixel.begin(), pixel.end(), kNaN);
        return pixel;
    }

    std::vector<double> ASDFgwcs::Convert(const std::string& from, const std::string& to, const std::vector<double>& coords) const
    {
        const size_t i = frameIndex(from);
        const size_t j = frameIndex(to);
        const size_t dims = i < fsteps.size() ? fsteps[i]->Inputs() : fsteps.back()->Outputs();
        if(coords.size() != dims)
            throw WCSexception(BAD_DIMEN, kWho, "Convert", "frame '"+from+"' has "+std::to_string(dims)+" axes, not "+std::to_string(coords.size()));
        return run(i, j, coords);
    }
}
