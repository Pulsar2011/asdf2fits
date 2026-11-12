//
//  ASDFtransform.cxx
//  asdf2fits
//
//  Written by Claude (Anthropic's AI model) for William Gillard.
//  Copyright (c) 2026 William Gillard.
//
//  SPDX-License-Identifier: BSD-3-Clause
//  Licence text: see the LICENSE file at the root of this project.

/**
 * @file ASDFtransform.cxx
 * @brief Implementation of DSL::ASDFtransform and of the kinds of transform known by the library.
 * @details The formulas are those of astropy.modeling and gwcs, which write these transforms; each class says
 * which model it follows. Angles are in degrees.
 */
#include "ASDFtransform.h"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <map>
#include <mutex>
#include <numbers>

#include <DSTfits/FITSimg.h>

#include "ASDFarray.h"

namespace DSL
{
    namespace
    {
        const char* const kWho = "ASDFtransform";   //!< Name given in the exceptions thrown here
        constexpr double kRad = std::numbers::pi / 180.;    //!< Radians per degree

        /**
         * @brief Room for the values passed between transforms, on the stack when there are few of them.
         */
        class Buffer
        {
        public:
            /** @brief Make room. @param size Number of values. */
            explicit Buffer(size_t size)
            {
                if(size > ffixed.size())
                    fheap.resize(size);
                fdata = size > ffixed.size() ? fheap.data() : ffixed.data();
            }
            double* data() { return fdata; }                        //!< @brief The values. @return Pointer to the first.
            double& operator[](size_t i) { return fdata[i]; }       //!< @brief A value. @param i Index. @return The value.

        private:
            std::array<double, 16> ffixed{};    //!< Room for few values
            std::vector<double> fheap;          //!< Room for many values
            double* fdata = nullptr;            //!< Where the values are
        };

        /**
         * @brief Map a value from a domain to a window, as astropy's poly_map_domain().
         * @param x      The value.
         * @param domain Domain [lower, upper].
         * @param window Window [lower, upper].
         * @return The value in the window.
         */
        double mapDomain(double x, const double* domain, const double* window)
        {
            const double scale  = (window[1] - window[0]) / (domain[1] - domain[0]);
            const double offset = (window[0]*domain[1] - window[1]*domain[0]) / (domain[1] - domain[0]);
            return offset + scale*x;
        }

        /**
         * @brief Read a whole number of a transform.
         * @param node     Node of the transform.
         * @param key      Key of the number.
         * @param fallback Value if the key is missing.
         * @return The number.
         */
        size_t count(const ASDFnode& node, const std::string& key, size_t fallback)
        {
            const double value = ASDFtransform::Number(node, key, static_cast<double>(fallback));
            if(value < 0 || value != std::floor(value))
                throw WCSexception(BAD_OPTION, kWho, "count", "'"+key+"' must be a whole number, not "+std::to_string(value));
            return static_cast<size_t>(value);
        }

        /**
         * @brief Read a text of a transform.
         * @param node     Node of the transform.
         * @param key      Key of the text.
         * @param fallback Value if the key is missing.
         * @return The text.
         */
        std::string text(const ASDFnode& node, const std::string& key, const std::string& fallback)
        {
            const ASDFnode* value = node.Find(key);
            if(!value || value->IsNull())
                return fallback;
            if(!value->IsString())
                throw WCSexception(BAD_OPTION, kWho, "text", "'"+key+"' must be a text");
            return value->Get<std::string>();
        }

#pragma region - Compositions
        /** @brief Second transform applied to the outputs of the first: astropy's `a | b`. */
        class Compose final : public ASDFtransform
        {
        public:
            /** @brief Chain two transforms. @param first Applied first. @param second Applied to its outputs. */
            Compose(std::shared_ptr<const ASDFtransform> first, std::shared_ptr<const ASDFtransform> second):
                ASDFtransform(first->Inputs(), second->Outputs()), ffirst(std::move(first)), fsecond(std::move(second))
            {
                if(ffirst->Outputs() != fsecond->Inputs())
                    Fail("compose", "the first transform gives "+std::to_string(ffirst->Outputs())+" values, the second takes "+std::to_string(fsecond->Inputs()));
            }

            void Evaluate(const double* in, double* out) const override
            {
                Buffer middle(ffirst->Outputs());
                ffirst->Evaluate(in, middle.data());
                fsecond->Evaluate(middle.data(), out);
            }

        protected:
            std::shared_ptr<const ASDFtransform> AnalyticInverse() const override
            {
                std::shared_ptr<const ASDFtransform> first = ffirst->Inverse();
                std::shared_ptr<const ASDFtransform> second = fsecond->Inverse();
                if(!first || !second)
                    return nullptr;
                return std::make_shared<Compose>(second, first);
            }

        private:
            std::shared_ptr<const ASDFtransform> ffirst;    //!< Applied first
            std::shared_ptr<const ASDFtransform> fsecond;   //!< Applied to the outputs of the first
        };

        /** @brief Transforms side by side, each on its own inputs: astropy's `a & b`. */
        class Concatenate final : public ASDFtransform
        {
        public:
            /** @brief Put transforms side by side. @param parts The transforms, first inputs first. */
            explicit Concatenate(std::vector<std::shared_ptr<const ASDFtransform>> parts):
                ASDFtransform(total(parts, true), total(parts, false)), fparts(std::move(parts)) {}

            void Evaluate(const double* in, double* out) const override
            {
                for(const auto& part : fparts)
                {
                    part->Evaluate(in, out);
                    in  += part->Inputs();
                    out += part->Outputs();
                }
            }

        protected:
            std::shared_ptr<const ASDFtransform> AnalyticInverse() const override
            {
                std::vector<std::shared_ptr<const ASDFtransform>> inverses;
                for(const auto& part : fparts)
                {
                    inverses.push_back(part->Inverse());
                    if(!inverses.back())
                        return nullptr;
                }
                return std::make_shared<Concatenate>(std::move(inverses));
            }

        private:
            /** @brief Add up inputs or outputs. @param parts The transforms. @param inputs Inputs or outputs. @return The total. */
            static size_t total(const std::vector<std::shared_ptr<const ASDFtransform>>& parts, bool inputs)
            {
                size_t sum = 0;
                for(const auto& part : parts)
                    sum += inputs ? part->Inputs() : part->Outputs();
                return sum;
            }

            std::vector<std::shared_ptr<const ASDFtransform>> fparts;   //!< The transforms, first inputs first
        };

        /** @brief Two transforms of the same inputs combined value by value: astropy's `a + b`, `a / b`... */
        class Arithmetic final : public ASDFtransform
        {
        public:
            enum class Op { add, subtract, multiply, divide, power };   //!< How values are combined

            /** @brief Combine two transforms. @param op How. @param left Left operand. @param right Right operand. */
            Arithmetic(Op op, std::shared_ptr<const ASDFtransform> left, std::shared_ptr<const ASDFtransform> right):
                ASDFtransform(left->Inputs(), left->Outputs()), fop(op), fleft(std::move(left)), fright(std::move(right))
            {
                if(fleft->Inputs() != fright->Inputs() || fleft->Outputs() != fright->Outputs())
                    Fail("arithmetic", "both operands must have the same inputs and outputs");
            }

            void Evaluate(const double* in, double* out) const override
            {
                Buffer right(Outputs());
                fleft->Evaluate(in, out);
                fright->Evaluate(in, right.data());
                for(size_t k = 0; k < Outputs(); k++)
                    switch(fop)
                    {
                        case Op::add:      out[k] += right[k]; break;
                        case Op::subtract: out[k] -= right[k]; break;
                        case Op::multiply: out[k] *= right[k]; break;
                        case Op::divide:   out[k] /= right[k]; break;
                        case Op::power:    out[k] = std::pow(out[k], right[k]); break;
                    }
            }

        private:
            Op fop;                                         //!< How values are combined
            std::shared_ptr<const ASDFtransform> fleft;     //!< Left operand
            std::shared_ptr<const ASDFtransform> fright;    //!< Right operand
        };

#pragma endregion
#pragma region - Simple transforms
        /** @brief Adds an offset: astropy's Shift. */
        class Shift final : public ASDFtransform
        {
        public:
            /** @brief Make the shift. @param offset Offset added. */
            explicit Shift(double offset): ASDFtransform(1, 1), foffset(offset) {}
            void Evaluate(const double* in, double* out) const override { out[0] = in[0] + foffset; }

        protected:
            std::shared_ptr<const ASDFtransform> AnalyticInverse() const override { return std::make_shared<Shift>(-foffset); }

        private:
            double foffset; //!< Offset added
        };

        /** @brief Multiplies by a factor: astropy's Scale and Multiply. */
        class Scale final : public ASDFtransform
        {
        public:
            /** @brief Make the scale. @param factor Factor. */
            explicit Scale(double factor): ASDFtransform(1, 1), ffactor(factor) {}
            void Evaluate(const double* in, double* out) const override { out[0] = in[0] * ffactor; }

        protected:
            std::shared_ptr<const ASDFtransform> AnalyticInverse() const override
            {
                return ffactor != 0. ? std::make_shared<Scale>(1./ffactor) : nullptr;
            }

        private:
            double ffactor; //!< Factor
        };

        /** @brief Gives its inputs back: astropy's Identity. */
        class Identity final : public ASDFtransform
        {
        public:
            /** @brief Make the identity. @param dims Number of inputs. */
            explicit Identity(size_t dims): ASDFtransform(dims, dims) {}
            void Evaluate(const double* in, double* out) const override { std::copy(in, in + Inputs(), out); }

        protected:
            std::shared_ptr<const ASDFtransform> AnalyticInverse() const override { return std::make_shared<Identity>(Inputs()); }
        };

        /** @brief Reorders, repeats or drops inputs: astropy's Mapping. */
        class Mapping final : public ASDFtransform
        {
        public:
            /**
             * @brief Make the mapping.
             * @param mapping Input given as each output.
             * @param inputs  Number of inputs.
             */
            Mapping(std::vector<size_t> mapping, size_t inputs): ASDFtransform(inputs, mapping.size()), fmapping(std::move(mapping))
            {
                for(size_t index : fmapping)
                    if(index >= inputs)
                        Fail("remap_axes", "input "+std::to_string(index)+" of "+std::to_string(inputs));
            }

            void Evaluate(const double* in, double* out) const override
            {
                for(size_t k = 0; k < fmapping.size(); k++)
                    out[k] = in[fmapping[k]];
            }

        protected:
            /** As astropy: each input is taken back from the first output it went to; none if one was dropped. */
            std::shared_ptr<const ASDFtransform> AnalyticInverse() const override
            {
                std::vector<size_t> back;
                for(size_t input = 0; input < Inputs(); input++)
                {
                    auto it = std::find(fmapping.begin(), fmapping.end(), input);
                    if(it == fmapping.end())
                        return nullptr;
                    back.push_back(static_cast<size_t>(it - fmapping.begin()));
                }
                return std::make_shared<Mapping>(std::move(back), fmapping.size());
            }

        private:
            std::vector<size_t> fmapping;   //!< Input given as each output
        };

        /** @brief Gives the same value whatever its inputs: astropy's Const1D and Const2D. */
        class Constant final : public ASDFtransform
        {
        public:
            /** @brief Make the constant. @param value Value given. @param inputs Number of inputs. */
            Constant(double value, size_t inputs): ASDFtransform(inputs, 1), fvalue(value) {}
            void Evaluate(const double*, double* out) const override { out[0] = fvalue; }

        private:
            double fvalue;  //!< Value given
        };

#pragma endregion
#pragma region - Polynomials and linear transforms
        /**
         * @brief Polynomial of one or two inputs: astropy's Polynomial1D and Polynomial2D.
         * @details In two inputs, the coefficient of x^i y^j is at [i, j]; only terms with i + j up to the degree
         * count. The inputs are first mapped from the domain to the window, when a domain is given.
         */
        class Polynomial final : public ASDFtransform
        {
        public:
            /**
             * @brief Make the polynomial.
             * @param inputs       1 or 2.
             * @param degree       Degree.
             * @param coefficients (degree+1) values for one input, (degree+1)^2 for two, C order.
             * @param domain       Domain of each input, [lower, upper]...; empty for none.
             * @param window       Window of each input, as @p domain.
             */
            Polynomial(size_t inputs, size_t degree, std::vector<double> coefficients, std::vector<double> domain, std::vector<double> window):
                ASDFtransform(inputs, 1), fdegree(degree), fcoefficients(std::move(coefficients)),
                fdomain(std::move(domain)), fwindow(std::move(window)) {}

            void Evaluate(const double* in, double* out) const override
            {
                double x = in[0];
                if(!fdomain.empty())
                    x = mapDomain(x, &fdomain[0], &fwindow[0]);

                if(Inputs() == 1)
                {
                    // Horner
                    double sum = 0.;
                    for(size_t i = fdegree + 1; i-- > 0;)
                        sum = sum*x + fcoefficients[i];
                    out[0] = sum;
                    return;
                }

                double y = in[1];
                if(!fdomain.empty())
                    y = mapDomain(y, &fdomain[2], &fwindow[2]);

                Buffer xPower(fdegree + 1), yPower(fdegree + 1);
                xPower[0] = yPower[0] = 1.;
                for(size_t i = 1; i <= fdegree; i++)
                {
                    xPower[i] = xPower[i-1]*x;
                    yPower[i] = yPower[i-1]*y;
                }

                const size_t n = fdegree + 1;
                double sum = 0.;
                for(size_t i = 0; i <= fdegree; i++)
                    for(size_t j = 0; i + j <= fdegree; j++)
                        sum += fcoefficients[i*n + j] * xPower[i] * yPower[j];
                out[0] = sum;
            }

        private:
            size_t fdegree;                     //!< Degree
            std::vector<double> fcoefficients;  //!< Coefficients, C order
            std::vector<double> fdomain;        //!< Domain of each input; empty for none
            std::vector<double> fwindow;        //!< Window of each input
        };

        /** @brief Matrix product plus translation: astropy's AffineTransformation2D, in any dimension. */
        class Affine final : public ASDFtransform
        {
        public:
            /** @brief Make the transform. @param matrix n x n values, C order. @param translation n values. */
            Affine(std::vector<double> matrix, std::vector<double> translation):
                ASDFtransform(translation.size(), translation.size()), fmatrix(std::move(matrix)), ftranslation(std::move(translation)) {}

            void Evaluate(const double* in, double* out) const override
            {
                const size_t n = ftranslation.size();
                for(size_t r = 0; r < n; r++)
                {
                    double sum = ftranslation[r];
                    for(size_t c = 0; c < n; c++)
                        sum += fmatrix[r*n + c] * in[c];
                    out[r] = sum;
                }
            }

        protected:
            /** The inverse matrix, by Gauss-Jordan elimination with partial pivoting; none if singular. */
            std::shared_ptr<const ASDFtransform> AnalyticInverse() const override
            {
                const size_t n = ftranslation.size();
                std::vector<double> a(fmatrix);
                std::vector<double> inv(n*n, 0.);
                for(size_t i = 0; i < n; i++)
                    inv[i*n + i] = 1.;

                for(size_t col = 0; col < n; col++)
                {
                    size_t pivot = col;
                    for(size_t r = col + 1; r < n; r++)
                        if(std::fabs(a[r*n + col]) > std::fabs(a[pivot*n + col]))
                            pivot = r;
                    if(a[pivot*n + col] == 0.)
                        return nullptr;
                    for(size_t c = 0; c < n; c++)
                    {
                        std::swap(a[col*n + c], a[pivot*n + c]);
                        std::swap(inv[col*n + c], inv[pivot*n + c]);
                    }
                    const double scale = a[col*n + col];
                    for(size_t c = 0; c < n; c++)
                    {
                        a[col*n + c] /= scale;
                        inv[col*n + c] /= scale;
                    }
                    for(size_t r = 0; r < n; r++)
                        if(r != col)
                        {
                            const double factor = a[r*n + col];
                            for(size_t c = 0; c < n; c++)
                            {
                                a[r*n + c] -= factor * a[col*n + c];
                                inv[r*n + c] -= factor * inv[col*n + c];
                            }
                        }
                }

                std::vector<double> back(n, 0.);
                for(size_t r = 0; r < n; r++)
                    for(size_t c = 0; c < n; c++)
                        back[r] -= inv[r*n + c] * ftranslation[c];
                return std::make_shared<Affine>(std::move(inv), std::move(back));
            }

        private:
            std::vector<double> fmatrix;        //!< n x n, C order
            std::vector<double> ftranslation;   //!< n values
        };

#pragma endregion
#pragma region - Spherical transforms
        /** @brief Unit vector of a direction. @param lon Longitude, degrees. @param lat Latitude, degrees. @param xyz Receives the vector. */
        void toCartesian(double lon, double lat, double* xyz)
        {
            const double c = std::cos(lat*kRad);
            xyz[0] = c*std::cos(lon*kRad);
            xyz[1] = c*std::sin(lon*kRad);
            xyz[2] = std::sin(lat*kRad);
        }

        /**
         * @brief Direction of a vector, as gwcs's CartesianToSpherical.
         * @param xyz    The vector.
         * @param lonLat Receives longitude and latitude, degrees.
         * @param wrap   180 for a longitude in (-180, 180], 360 for [0, 360).
         */
        void toSpherical(const double* xyz, double* lonLat, double wrap)
        {
            lonLat[1] = std::atan2(xyz[2], std::hypot(xyz[0], xyz[1])) / kRad;
            double lon = std::atan2(xyz[1], xyz[0]) / kRad;
            if(wrap != 180.)
            {
                lon = std::fmod(lon, 360.);
                if(lon < 0.)
                    lon += 360.;
            }
            lonLat[0] = lon;
        }

        /** @brief Between longitude and latitude and a unit vector: gwcs's SphericalToCartesian and CartesianToSpherical. */
        class SphericalCartesian final : public ASDFtransform
        {
        public:
            /**
             * @brief Make the transform.
             * @param toCartesianWay From longitude and latitude to a vector, or back.
             * @param wrap          Wrapping of the longitude given back: 180 or 360.
             */
            SphericalCartesian(bool toCartesianWay, double wrap):
                ASDFtransform(toCartesianWay ? 2 : 3, toCartesianWay ? 3 : 2), fcartesian(toCartesianWay), fwrap(wrap) {}

            void Evaluate(const double* in, double* out) const override
            {
                if(fcartesian)
                    toCartesian(in[0], in[1], out);
                else
                    toSpherical(in, out, fwrap);
            }

        protected:
            std::shared_ptr<const ASDFtransform> AnalyticInverse() const override
            {
                return std::make_shared<SphericalCartesian>(!fcartesian, fwrap);
            }

        private:
            bool fcartesian;    //!< From longitude and latitude to a vector
            double fwrap;       //!< Wrapping of the longitude given back
        };

        /**
         * @brief Rotations about the axes, one after the other: astropy's RotationSequence3D and SphericalRotationSequence.
         * @details Each rotation turns the frame (astropy's rotation_matrix()), the first angle first. The cartesian
         * form works on vectors, the spherical one on longitude and latitude.
         */
        class RotationSequence final : public ASDFtransform
        {
        public:
            /**
             * @brief Make the rotation.
             * @param angles    Angles, degrees.
             * @param axes      Axis of each angle: 'x', 'y' or 'z'.
             * @param spherical Work on longitude and latitude rather than vectors.
             */
            RotationSequence(std::vector<double> angles, std::string axes, bool spherical):
                ASDFtransform(spherical ? 2 : 3, spherical ? 2 : 3), fangles(std::move(angles)), faxes(std::move(axes)), fspherical(spherical)
            {
                if(fangles.size() != faxes.size())
                    Fail("rotate_sequence_3d", std::to_string(fangles.size())+" angles for "+std::to_string(faxes.size())+" axes");

                fmatrix = {1., 0., 0., 0., 1., 0., 0., 0., 1.};
                for(size_t k = 0; k < fangles.size(); k++)
                {
                    const char axis = static_cast<char>(std::tolower(static_cast<unsigned char>(faxes[k])));
                    if(axis < 'x' || axis > 'z')
                        Fail("rotate_sequence_3d", std::string("unknown axis '")+faxes[k]+"'");
                    const size_t i = static_cast<size_t>(axis - 'x');
                    const size_t a1 = (i + 1) % 3, a2 = (i + 2) % 3;
                    const double c = std::cos(fangles[k]*kRad), s = std::sin(fangles[k]*kRad);

                    std::array<double, 9> r{};
                    r[i*3 + i] = 1.;
                    r[a1*3 + a1] = c;
                    r[a1*3 + a2] = s;
                    r[a2*3 + a1] = -s;
                    r[a2*3 + a2] = c;

                    // This rotation applies after the previous ones: matrix = r . matrix
                    std::array<double, 9> product{};
                    for(size_t row = 0; row < 3; row++)
                        for(size_t col = 0; col < 3; col++)
                            for(size_t m = 0; m < 3; m++)
                                product[row*3 + col] += r[row*3 + m] * fmatrix[m*3 + col];
                    fmatrix = product;
                }
            }

            void Evaluate(const double* in, double* out) const override
            {
                double xyz[3];
                if(fspherical)
                    toCartesian(in[0], in[1], xyz);
                else
                    std::copy(in, in + 3, xyz);

                double turned[3];
                for(size_t row = 0; row < 3; row++)
                    turned[row] = fmatrix[row*3]*xyz[0] + fmatrix[row*3 + 1]*xyz[1] + fmatrix[row*3 + 2]*xyz[2];

                if(fspherical)
                    toSpherical(turned, out, 180.);
                else
                    std::copy(turned, turned + 3, out);
            }

        protected:
            /** As astropy: the opposite angles, in reverse order. */
            std::shared_ptr<const ASDFtransform> AnalyticInverse() const override
            {
                std::vector<double> angles(fangles.rbegin(), fangles.rend());
                for(double& angle : angles)
                    angle = -angle;
                return std::make_shared<RotationSequence>(std::move(angles), std::string(faxes.rbegin(), faxes.rend()), fspherical);
            }

        private:
            std::vector<double> fangles;        //!< Angles, degrees
            std::string faxes;                  //!< Axis of each angle
            bool fspherical;                    //!< Works on longitude and latitude
            std::array<double, 9> fmatrix{};    //!< Product of the rotations, row by row
        };

        /**
         * @brief Gnomonic (TAN) projection: astropy's Pix2Sky_Gnomonic and Sky2Pix_Gnomonic.
         * @details Between projection plane coordinates (x, y), degrees, and native spherical coordinates (phi, theta).
         */
        class Gnomonic final : public ASDFtransform
        {
        public:
            /** @brief Make the projection. @param toSky From the plane to the sphere, or back. */
            explicit Gnomonic(bool toSky): ASDFtransform(2, 2), fsky(toSky) {}

            void Evaluate(const double* in, double* out) const override
            {
                if(fsky)
                {
                    const double r = std::hypot(in[0], in[1]);
                    out[0] = r == 0. ? 0. : std::atan2(in[0], -in[1]) / kRad;
                    out[1] = std::atan2(1./kRad, r) / kRad;
                }
                else
                {
                    const double r = std::cos(in[1]*kRad) / std::sin(in[1]*kRad) / kRad;
                    out[0] =  r*std::sin(in[0]*kRad);
                    out[1] = -r*std::cos(in[0]*kRad);
                }
            }

        protected:
            std::shared_ptr<const ASDFtransform> AnalyticInverse() const override { return std::make_shared<Gnomonic>(!fsky); }

        private:
            bool fsky;  //!< From the plane to the sphere
        };

        /**
         * @brief Rotation between native and celestial spherical coordinates: astropy's RotateNative2Celestial
         * and RotateCelestial2Native, as in Calabretta & Greisen (2002), equations 2 and 5.
         * @details The celestial longitude given back is in [0, 360).
         */
        class SkyRotation final : public ASDFtransform
        {
        public:
            /**
             * @brief Make the rotation.
             * @param lon         Celestial longitude of the native pole (alpha_p), degrees.
             * @param lat         Celestial latitude of the native pole (delta_p), degrees.
             * @param lonPole     Native longitude of the celestial pole (phi_p), degrees.
             * @param toCelestial From native to celestial coordinates, or back.
             */
            SkyRotation(double lon, double lat, double lonPole, bool toCelestial):
                ASDFtransform(2, 2), flon(lon), flat(lat), flonPole(lonPole), fcelestial(toCelestial) {}

            void Evaluate(const double* in, double* out) const override
            {
                const double sinP = std::sin(flat*kRad), cosP = std::cos(flat*kRad);
                const double sinB = std::sin(in[1]*kRad), cosB = std::cos(in[1]*kRad);
                // From one pole to the other: native (phi, theta) <-> celestial (alpha, delta) share these formulas
                const double d = (in[0] - (fcelestial ? flonPole : flon)) * kRad;
                double lon = (fcelestial ? flon : flonPole)
                           + std::atan2(-cosB*std::sin(d), sinB*cosP - cosB*sinP*std::cos(d)) / kRad;
                const double lat = std::asin(std::clamp(sinB*sinP + cosB*cosP*std::cos(d), -1., 1.)) / kRad;
                lon = std::fmod(lon, 360.);
                if(lon < 0.)
                    lon += 360.;
                out[0] = lon;
                out[1] = lat;
            }

        protected:
            std::shared_ptr<const ASDFtransform> AnalyticInverse() const override
            {
                return std::make_shared<SkyRotation>(flon, flat, flonPole, !fcelestial);
            }

        private:
            double flon;        //!< Celestial longitude of the native pole
            double flat;        //!< Celestial latitude of the native pole
            double flonPole;    //!< Native longitude of the celestial pole
            bool fcelestial;    //!< From native to celestial coordinates
        };

#pragma endregion
#pragma region - Builders
        /** @brief Read the two operands of a composition or an arithmetic operator. @param node Its node. @param where Kind. @return The operands. */
        std::vector<std::shared_ptr<const ASDFtransform>> twoOperands(const ASDFnode& node, const std::string& where)
        {
            std::vector<std::shared_ptr<const ASDFtransform>> parts = ASDFtransform::Operands(node);
            if(parts.size() < 2)
                throw WCSexception(BAD_OPTION, kWho, where, "needs two transforms under 'forward', not "+std::to_string(parts.size()));
            return parts;
        }

        /** @brief Register the kinds of transform known by the library. @param builders Where to register them. */
        void addBuiltins(std::map<std::string, ASDFtransform::Builder>& builders)
        {
            builders["transform/compose"] = [](const ASDFnode& node) -> std::shared_ptr<ASDFtransform>
            {
                const auto parts = twoOperands(node, "compose");
                std::shared_ptr<const ASDFtransform> chain = parts[0];
                for(size_t k = 1; k + 1 < parts.size(); k++)
                    chain = std::make_shared<Compose>(chain, parts[k]);
                return std::make_shared<Compose>(chain, parts.back());
            };

            builders["transform/concatenate"] = [](const ASDFnode& node) -> std::shared_ptr<ASDFtransform>
            {
                return std::make_shared<Concatenate>(ASDFtransform::Operands(node));
            };

            const std::pair<const char*, Arithmetic::Op> operators[] = {
                {"transform/add", Arithmetic::Op::add},       {"transform/subtract", Arithmetic::Op::subtract},
                {"transform/multiply", Arithmetic::Op::multiply}, {"transform/divide", Arithmetic::Op::divide},
                {"transform/power", Arithmetic::Op::power}};
            for(const auto& entry : operators)
            {
                const std::string name = entry.first;
                const Arithmetic::Op op = entry.second;
                builders[name] = [op, name](const ASDFnode& node) -> std::shared_ptr<ASDFtransform>
                {
                    const auto parts = twoOperands(node, name);
                    return std::make_shared<Arithmetic>(op, parts[0], parts[1]);
                };
            }

            builders["transform/shift"] = [](const ASDFnode& node) -> std::shared_ptr<ASDFtransform>
            {
                return std::make_shared<Shift>(ASDFtransform::Number(node, "offset"));
            };

            builders["transform/scale"] = [](const ASDFnode& node) -> std::shared_ptr<ASDFtransform>
            {
                return std::make_shared<Scale>(ASDFtransform::Number(node, "factor"));
            };
            builders["transform/multiplyscale"] = builders["transform/scale"];

            builders["transform/identity"] = [](const ASDFnode& node) -> std::shared_ptr<ASDFtransform>
            {
                return std::make_shared<Identity>(count(node, "n_dims", 1));
            };

            builders["transform/remap_axes"] = [](const ASDFnode& node) -> std::shared_ptr<ASDFtransform>
            {
                std::vector<size_t> mapping;
                size_t inputs = 0;
                for(double index : ASDFtransform::Numbers(ASDFtransform::Required(node, "mapping")))
                {
                    if(index < 0 || index != std::floor(index))
                        throw WCSexception(BAD_OPTION, kWho, "remap_axes", "the mapping must hold input indices");
                    mapping.push_back(static_cast<size_t>(index));
                    inputs = std::max(inputs, mapping.back() + 1);
                }
                return std::make_shared<Mapping>(std::move(mapping), count(node, "n_inputs", inputs));
            };

            builders["transform/constant"] = [](const ASDFnode& node) -> std::shared_ptr<ASDFtransform>
            {
                return std::make_shared<Constant>(ASDFtransform::Number(node, "value"), count(node, "dimensions", 1));
            };

            builders["transform/polynomial"] = [](const ASDFnode& node) -> std::shared_ptr<ASDFtransform>
            {
                const ASDFnode& coefficients = ASDFtransform::Required(node, "coefficients");
                const std::vector<size_t> shape = ASDFtransform::ShapeOf(coefficients);
                std::vector<double> values = ASDFtransform::Numbers(coefficients);
                if(shape.empty() || shape.size() > 2 || (shape.size() == 2 && shape[0] != shape[1]) || shape[0] == 0)
                    throw WCSexception(BAD_OPTION, kWho, "polynomial", "the coefficients must be a vector or a square matrix");
                const size_t inputs = shape.size();

                std::vector<double> domain, window;
                if(const ASDFnode* d = node.Find("domain"); d && !d->IsNull())
                {
                    domain = ASDFtransform::Numbers(*d);
                    const ASDFnode* w = node.Find("window");
                    window = (w && !w->IsNull()) ? ASDFtransform::Numbers(*w) : std::vector<double>();
                    if(window.empty())
                        for(size_t k = 0; k < inputs; k++)
                            window.insert(window.end(), {-1., 1.});
                    if(domain.size() != 2*inputs || window.size() != 2*inputs)
                        throw WCSexception(BAD_OPTION, kWho, "polynomial", "domain and window need two values per input");
                    for(size_t k = 0; k < inputs; k++)
                        if(domain[2*k] == domain[2*k + 1])
                            throw WCSexception(BAD_OPTION, kWho, "polynomial", "empty domain");
                }
                return std::make_shared<Polynomial>(inputs, shape[0] - 1, std::move(values), std::move(domain), std::move(window));
            };

            builders["transform/affine"] = [](const ASDFnode& node) -> std::shared_ptr<ASDFtransform>
            {
                const ASDFnode& matrix = ASDFtransform::Required(node, "matrix");
                std::vector<double> m = ASDFtransform::Numbers(matrix);
                const std::vector<size_t> shape = ASDFtransform::ShapeOf(matrix);
                std::vector<double> t;
                if(const ASDFnode* translation = node.Find("translation"); translation && !translation->IsNull())
                    t = ASDFtransform::Numbers(*translation);
                else if(shape.size() == 2)
                    t.assign(shape[0], 0.);
                if(shape.size() != 2 || shape[0] != shape[1] || t.size() != shape[0])
                    throw WCSexception(BAD_OPTION, kWho, "affine", "the matrix must be square and the translation as long as a side");
                return std::make_shared<Affine>(std::move(m), std::move(t));
            };

            builders["transform/rotate_sequence_3d"] = [](const ASDFnode& node) -> std::shared_ptr<ASDFtransform>
            {
                const std::string type = text(node, "rotation_type", "cartesian");
                if(type != "cartesian" && type != "spherical")
                    throw WCSexception(BAD_OPTION, kWho, "rotate_sequence_3d", "unknown rotation_type '"+type+"'");
                return std::make_shared<RotationSequence>(ASDFtransform::Numbers(ASDFtransform::Required(node, "angles")),
                                                          text(node, "axes_order", ""), type == "spherical");
            };

            builders["gwcs/spherical_cartesian"] = [](const ASDFnode& node) -> std::shared_ptr<ASDFtransform>
            {
                const std::string type = text(node, "transform_type", "");
                if(type != "spherical_to_cartesian" && type != "cartesian_to_spherical")
                    throw WCSexception(BAD_OPTION, kWho, "spherical_cartesian", "unknown transform_type '"+type+"'");
                const double wrap = ASDFtransform::Number(node, "wrap_lon_at", 360.);
                if(wrap != 180. && wrap != 360.)
                    throw WCSexception(BAD_OPTION, kWho, "spherical_cartesian", "wrap_lon_at must be 180 or 360");
                return std::make_shared<SphericalCartesian>(type == "spherical_to_cartesian", wrap);
            };

            builders["transform/gnomonic"] = [](const ASDFnode& node) -> std::shared_ptr<ASDFtransform>
            {
                const std::string direction = text(node, "direction", "pix2sky");
                if(direction != "pix2sky" && direction != "sky2pix")
                    throw WCSexception(BAD_OPTION, kWho, "gnomonic", "unknown direction '"+direction+"'");
                return std::make_shared<Gnomonic>(direction == "pix2sky");
            };

            builders["transform/rotate3d"] = [](const ASDFnode& node) -> std::shared_ptr<ASDFtransform>
            {
                const std::string direction = text(node, "direction", "native2celestial");
                if(direction != "native2celestial" && direction != "celestial2native")
                    throw WCSexception(BAD_OPTION, kWho, "rotate3d", "direction '"+direction+"' (Euler angles) isn't supported");
                return std::make_shared<SkyRotation>(ASDFtransform::Number(node, "phi"), ASDFtransform::Number(node, "theta"),
                                                     ASDFtransform::Number(node, "psi"), direction == "native2celestial");
            };
        }

        /** @brief The builders, by kind. */
        struct Registry
        {
            Registry() { addBuiltins(builders); }   //!< Starts with the kinds known by the library
            std::mutex lock;                                        //!< Guards the builders
            std::map<std::string, ASDFtransform::Builder> builders; //!< Builder of each kind
        };

        /** @brief Get the builders. @return The registry, made on first use. */
        Registry& registry()
        {
            static Registry instance;
            return instance;
        }

#pragma endregion
    }

#pragma region - ASDFtransform
    ASDFtransform::ASDFtransform(size_t inputs, size_t outputs): finputs(inputs), foutputs(outputs) {}

    ASDFtransform::~ASDFtransform() = default;

    void ASDFtransform::Fail(const std::string& where, const std::string& what)
    {
        throw WCSexception(BAD_OPTION, kWho, where, what);
    }

    std::vector<double> ASDFtransform::operator()(const std::vector<double>& in) const
    {
        if(in.size() != finputs)
            throw WCSexception(BAD_DIMEN, kWho, "operator()", "the transform takes "+std::to_string(finputs)+" values, not "+std::to_string(in.size()));
        std::vector<double> out(foutputs);
        Evaluate(in.data(), out.data());
        return out;
    }

    std::shared_ptr<const ASDFtransform> ASDFtransform::Inverse() const
    {
        return finverse ? finverse : AnalyticInverse();
    }

    std::string ASDFtransform::KindOf(const std::string& tag)
    {
        std::string name = tag;
        const size_t dash = name.rfind('-');
        if(dash != std::string::npos && dash + 1 < name.size() && std::isdigit(static_cast<unsigned char>(name[dash + 1])))
            name.resize(dash);

        std::vector<std::string> parts;
        std::string part;
        for(char c : name + "/")
            if(c == '/' || c == ':')
            {
                if(!part.empty() && part != "tags")
                    parts.push_back(part);
                part.clear();
            }
            else
                part += c;

        if(parts.empty())
            return std::string();
        if(parts.size() == 1)
            return parts[0];
        return parts[parts.size() - 2] + "/" + parts.back();
    }

    void ASDFtransform::Register(const std::string& kind, Builder builder)
    {
        Registry& r = registry();
        std::lock_guard<std::mutex> guard(r.lock);
        r.builders[kind] = std::move(builder);
    }

    bool ASDFtransform::IsKnown(const std::string& kind)
    {
        Registry& r = registry();
        std::lock_guard<std::mutex> guard(r.lock);
        return r.builders.count(kind) > 0;
    }

    std::vector<std::string> ASDFtransform::KnownKinds()
    {
        Registry& r = registry();
        std::lock_guard<std::mutex> guard(r.lock);
        std::vector<std::string> kinds;
        for(const auto& entry : r.builders)
            kinds.push_back(entry.first);
        return kinds;
    }

    std::shared_ptr<const ASDFtransform> ASDFtransform::Read(const ASDFnode& node)
    {
        if(!node.IsMapping())
            Fail("Read", "a transform must be a mapping");
        const std::string kind = KindOf(node.Tag());
        if(kind.empty())
            Fail("Read", "a transform must have a tag");

        Builder builder;
        {
            Registry& r = registry();
            std::lock_guard<std::mutex> guard(r.lock);
            auto it = r.builders.find(kind);
            if(it != r.builders.end())
                builder = it->second;
        }
        if(!builder)
            Fail("Read", "unknown kind of transform '"+kind+"' (tag "+node.Tag()+"); see ASDFtransform::Register()");

        std::shared_ptr<ASDFtransform> made = builder(node);
        if(!made)
            Fail("Read", "the builder of '"+kind+"' made nothing");
        made->fkind = kind;
        if(const ASDFnode* name = node.Find("name"); name && name->IsString())
            made->fname = name->Get<std::string>();

        // An inverse whose sizes don't fit is ignored, not an error: astropy writes some that it never evaluates,
        // such as identity(2) for the mapping (0, 1, 0, 1) in the JWST MIRI imager, inside a compose that has
        // its own inverse. The analytic inverse, or Newton's method in ASDFgwcs, is used instead.
        if(const ASDFnode* inverse = node.Find("inverse"); inverse && !inverse->IsNull())
        {
            std::shared_ptr<const ASDFtransform> back = Read(*inverse);
            if(back->Inputs() == made->foutputs && back->Outputs() == made->finputs)
                made->finverse = std::move(back);
        }
        return made;
    }

#pragma endregion
#pragma region - Reading of parameters
    const ASDFnode& ASDFtransform::Required(const ASDFnode& node, const std::string& key)
    {
        const ASDFnode* value = node.Find(key);
        if(!value || value->IsNull())
            Fail("Required", "'"+key+"' is missing from "+KindOf(node.Tag()));
        return *value;
    }

    double ASDFtransform::Number(const ASDFnode& node, const std::string& key)
    {
        const std::vector<double> values = Numbers(Required(node, key));
        if(values.size() != 1)
            Fail("Number", "'"+key+"' must be one number");
        return values[0];
    }

    double ASDFtransform::Number(const ASDFnode& node, const std::string& key, double fallback)
    {
        const ASDFnode* value = node.Find(key);
        if(!value || value->IsNull())
            return fallback;
        return Number(node, key);
    }

    std::vector<double> ASDFtransform::Numbers(const ASDFnode& node)
    {
        if(node.IsNumber())
            return {node.Get<double>()};

        if(node.IsArray())
        {
            // The whole array in one read; its elements come in C order, the last axis varying fastest
            std::shared_ptr<FITScube> cube = node.Array()->Read({}, {});
            std::vector<double> values(static_cast<size_t>(node.Array()->Elements()));
            for(size_t k = 0; k < values.size(); k++)
                values[k] = cube->at(k);
            return values;
        }

        if(node.IsSequence())
        {
            std::vector<double> values;
            for(const ASDFnode& item : node.Items())
            {
                const std::vector<double> more = Numbers(item);
                values.insert(values.end(), more.begin(), more.end());
            }
            return values;
        }

        if(node.IsMapping())
        {
            // A quantity, or an ndarray written inline
            if(const ASDFnode* value = node.Find("value"))
                return Numbers(*value);
            if(const ASDFnode* data = node.Find("data"))
                return Numbers(*data);
        }
        Fail("Numbers", "expected numbers, found "+(node.IsString() ? "'"+node.Text()+"'" : std::string("another kind of value")));
    }

    std::vector<size_t> ASDFtransform::ShapeOf(const ASDFnode& node)
    {
        if(node.IsArray())
        {
            const std::vector<uint64_t> shape = node.Array()->Shape();
            return std::vector<size_t>(shape.begin(), shape.end());
        }
        if(node.IsSequence())
        {
            std::vector<size_t> shape{node.Size()};
            if(node.Size() > 0)
            {
                const std::vector<size_t> inner = ShapeOf(node.Items()[0]);
                shape.insert(shape.end(), inner.begin(), inner.end());
            }
            return shape;
        }
        if(node.IsMapping())
        {
            if(const ASDFnode* value = node.Find("value"))
                return ShapeOf(*value);
            if(const ASDFnode* data = node.Find("data"))
                return ShapeOf(*data);
        }
        return {};
    }

    std::vector<std::shared_ptr<const ASDFtransform>> ASDFtransform::Operands(const ASDFnode& node)
    {
        const ASDFnode& forward = Required(node, "forward");
        if(!forward.IsSequence() || forward.Size() == 0)
            Fail("Operands", "'forward' of "+KindOf(node.Tag())+" must be a list of transforms");
        std::vector<std::shared_ptr<const ASDFtransform>> parts;
        for(const ASDFnode& item : forward.Items())
            parts.push_back(Read(item));
        return parts;
    }
#pragma endregion
}
