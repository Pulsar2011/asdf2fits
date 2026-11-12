//
//  ASDFtransform.h
//  asdf2fits
//
//  Written by Claude (Anthropic's AI model) for William Gillard.
//  Copyright (c) 2026 William Gillard.
//
//  SPDX-License-Identifier: BSD-3-Clause
//  Licence text: see the LICENSE file at the root of this project.

/**
 * @file ASDFtransform.h
 * @brief Transforms of a GWCS: DSL::ASDFtransform, the base class every kind of transform derives from.
 */
#ifndef _DSL_ASDFtransform_
#define _DSL_ASDFtransform_

#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "ASDFnode.h"

namespace DSL
{
    /**
     * @class DSL::ASDFtransform ASDFtransform.h "asdf2fits/ASDFtransform.h"
     * @brief One transform of a GWCS: a function from Inputs() numbers to Outputs() numbers.
     * @details A GWCS is a chain of transforms read from the tree of an ASDF file, where each transform is a
     * mapping tagged with its kind: "tag:stsci.edu:asdf/transform/shift-1.3.0", "tag:stsci.edu:gwcs/spherical_cartesian-1.2.0"...
     * Each kind is a class derived from ASDFtransform, made by a builder registered under the name of the kind
     * (see Register()). Read() finds the builder from the tag, so a kind of transform not yet known is supported
     * by deriving a class and registering its builder, from the library or from the code using it:
     *
     * @code
     * class Square : public DSL::ASDFtransform
     * {
     * public:
     *     Square(): ASDFtransform(1, 1) {}
     *     void Evaluate(const double* in, double* out) const override { out[0] = in[0]*in[0]; }
     * };
     *
     * DSL::ASDFtransform::Register("transform/square",
     *     [](const DSL::ASDFnode&) { return std::make_shared<Square>(); });
     * @endcode
     *
     * A transform can have an inverse. The file can give one, under the key "inverse" of the transform: Read() then
     * attaches it, and it comes first (unless its numbers of inputs and outputs don't fit, as some astropy writes and
     * never uses: it is then ignored). Otherwise Inverse() asks the transform for its own (AnalyticInverse()),
     * which a derived class gives when it has one: a shift gives the opposite shift, a composition the composition
     * of the inverses in reverse order. Without either, Inverse() is null and the caller has to invert numerically.
     *
     * Transforms are immutable once read, so they can be shared and used from several threads.
     * Angles are in degrees, as in astropy.modeling. Failures throw WCSexception.
     */
    class ASDFtransform
    {
    public:
        /**
         * @brief Make a transform from the node of a tree.
         * @details Takes the builder registered under the name of the tag of @p node (see KindOf()), and attaches
         * the inverse given under the key "inverse", if any. The key "name" names the transform.
         * @param node A mapping tagged as a transform.
         * @return The transform.
         * @throws WCSexception BAD_OPTION if no builder is registered for this kind, or as the builder.
         */
        static std::shared_ptr<const ASDFtransform> Read(const ASDFnode& node);

        /**
         * @brief Make a transform from a node of a tree.
         * @param node The node of the transform.
         * @return The new transform; Read() attaches its name and its inverse.
         * @throws WCSexception BAD_OPTION if the node doesn't describe a valid transform of this kind.
         */
        using Builder = std::function<std::shared_ptr<ASDFtransform>(const ASDFnode& node)>;

        /**
         * @brief Make a kind of transform known to Read().
         * @details Replaces the builder registered under the same name, if any: a builder of the library can be
         * overridden.
         * @param kind    Name of the kind, as KindOf() gives it: "transform/shift", "gwcs/spherical_cartesian"...
         * @param builder Makes a transform of this kind from its node.
         */
        static void Register(const std::string& kind, Builder builder);

        /**
         * @brief Tell whether a kind of transform is known.
         * @param kind Name of the kind, as KindOf() gives it.
         * @return True if a builder is registered under this name.
         */
        static bool IsKnown(const std::string& kind);

        /**
         * @brief Get the names of the kinds of transform known.
         * @return The names, in alphabetical order.
         */
        static std::vector<std::string> KnownKinds();

        /**
         * @brief Get the kind of a transform from its tag.
         * @details The kind is the name of the tag without its version, after the family it belongs to:
         * "tag:stsci.edu:asdf/transform/shift-1.3.0" and "asdf://asdf-format.org/transform/tags/shift-1.4.0" give
         * "transform/shift", "tag:stsci.edu:gwcs/spherical_cartesian-1.2.0" gives "gwcs/spherical_cartesian".
         * @param tag Tag of the node, in full.
         * @return The kind; empty if the tag is empty.
         */
        static std::string KindOf(const std::string& tag);

        virtual ~ASDFtransform();

        /** @brief Get the number of inputs. @return The number of values Evaluate() reads. */
        size_t Inputs() const { return finputs; }

        /** @brief Get the number of outputs. @return The number of values Evaluate() writes. */
        size_t Outputs() const { return foutputs; }

        /** @brief Get the name given to the transform in the file. @return The name; empty if none. */
        const std::string& Name() const { return fname; }

        /** @brief Get the kind of the transform. @return The kind it was read as, such as "transform/shift"; empty if made in code. */
        const std::string& Kind() const { return fkind; }

        /**
         * @brief Compute the outputs of the transform for one point.
         * @param in  Inputs() values.
         * @param out Receives Outputs() values; never the same memory as @p in.
         */
        virtual void Evaluate(const double* in, double* out) const = 0;

        /**
         * @brief Compute the outputs of the transform for one point.
         * @param in Inputs() values.
         * @return Outputs() values.
         * @throws WCSexception BAD_DIMEN if @p in hasn't Inputs() values.
         */
        std::vector<double> operator()(const std::vector<double>& in) const;

        /**
         * @brief Get the inverse.
         * @details The one given by the file under the key "inverse" if any; otherwise the one the transform
         * gives itself (AnalyticInverse()).
         * @return The inverse; null if there is none.
         */
        std::shared_ptr<const ASDFtransform> Inverse() const;

#pragma region * Reading of parameters, for the builders
        /**
         * @brief Read a number of a transform.
         * @details A quantity (a mapping with a "value") gives its value; its unit is ignored.
         * @param node Node of the transform.
         * @param key  Key of the number.
         * @return The number.
         * @throws WCSexception BAD_OPTION if the key is missing or isn't a number.
         */
        static double Number(const ASDFnode& node, const std::string& key);

        /**
         * @brief Read a number of a transform, if it is there.
         * @param node     Node of the transform.
         * @param key      Key of the number.
         * @param fallback Value if the key is missing.
         * @return The number, or @p fallback.
         * @throws WCSexception BAD_OPTION if the key is there but isn't a number.
         */
        static double Number(const ASDFnode& node, const std::string& key, double fallback);

        /**
         * @brief Read numbers: an ndarray, a sequence (of sequences...), a quantity or a single number.
         * @details An ndarray is read from its file, in the order of its elements (the last axis varies fastest).
         * @param node The node.
         * @return The numbers.
         * @throws WCSexception BAD_OPTION if the node holds something else.
         */
        static std::vector<double> Numbers(const ASDFnode& node);

        /**
         * @brief Get the shape of numbers read by Numbers().
         * @param node The node.
         * @return The shape, the slowest varying axis first: that of an ndarray, or the lengths of nested
         * sequences; empty for a single number.
         */
        static std::vector<size_t> ShapeOf(const ASDFnode& node);

        /**
         * @brief Get a value of a transform that must be there.
         * @param node Node of the transform.
         * @param key  Key of the value.
         * @return The value.
         * @throws WCSexception BAD_OPTION if the key is missing.
         */
        static const ASDFnode& Required(const ASDFnode& node, const std::string& key);

        /**
         * @brief Read the transforms listed under "forward", as compositions and concatenations have them.
         * @param node Node of the transform.
         * @return The transforms, in order.
         * @throws WCSexception BAD_OPTION if "forward" isn't a list of transforms.
         */
        static std::vector<std::shared_ptr<const ASDFtransform>> Operands(const ASDFnode& node);

#pragma endregion
    protected:
        /**
         * @brief Make a transform.
         * @param inputs  Number of inputs.
         * @param outputs Number of outputs.
         */
        ASDFtransform(size_t inputs, size_t outputs);

        /**
         * @brief Give the inverse a transform knows by itself.
         * @details Called by Inverse() when the file gives none. The default gives none.
         * @return The inverse; null if there is none.
         */
        virtual std::shared_ptr<const ASDFtransform> AnalyticInverse() const { return nullptr; }

        /**
         * @brief Throw the exception of a transform that can't be read.
         * @param where Builder or function failing.
         * @param what  What is wrong.
         * @throws WCSexception BAD_OPTION, always.
         */
        [[noreturn]] static void Fail(const std::string& where, const std::string& what);

    private:
        size_t finputs;                                 //!< Number of inputs
        size_t foutputs;                                //!< Number of outputs
        std::string fname;                              //!< Name given in the file
        std::string fkind;                              //!< Kind it was read as
        std::shared_ptr<const ASDFtransform> finverse;  //!< Inverse given by the file
    };
}

#endif
