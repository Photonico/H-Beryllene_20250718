#include "SphericalHarmonics.hpp"

#include "ParallelEnvironemt.hpp"
#include "constants.hpp"
#include "debug.hpp"
#include "math.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>

using namespace vaspml;
using Integrate3SphericalHarmonicsData = math::SphericalHarmonics::Integrate3SphericalHarmonicsData;

math::SphericalHarmonics::SphericalHarmonics(Int lmax) : ldim((lmax + 1) * (lmax + 1))
{
    integrate3SphericalHarmonics(lmax, this->i3sh);
}

VASPML_NV_HOST_DEVICE
void math::SphericalHarmonics::computeSphericalHarmonicsAndGradientsNormalized(const Vec1Real& rhat,
                                                                               Vec1Real&       ylm,
                                                                               Vec1Real& ylmd) const
{
    if (this->i3sh.lmax < 0) return;
    const Real        pre = 1.0 / std::sqrt(constants::PI4);
    const std::size_t n = rhat.size() / 3;

    // l = 0
    for (std::size_t i = 0; i < n * ldim; i += ldim)
    {
        ylm[i] = pre;
        ylmd[3 * i] = 0.0;
        ylmd[3 * i + 1] = 0.0;
        ylmd[3 * i + 2] = 0.0;
    }
    if (this->i3sh.lmax < 1) return;

    // l = 1
    const Real pres3 = std::sqrt(3.0 / constants::PI4);
    for (std::size_t i = 0; i < n; ++i)
    {
        std::size_t       s = i * ldim + 1;
        const std::size_t a = 3 * i;
        // clang-format off
        ylm[s    ] = pres3 * rhat[a + 1]; // l = 1, m = -1
        ylm[s + 1] = pres3 * rhat[a + 2]; // l = 1, m = 0
        ylm[s + 2] = pres3 * rhat[a    ]; // l = 1, m = 1
        s = 3 * i * ldim + 3;
        // l = 1, m = -1
        ylmd[s    ] = 0.0;
        ylmd[s + 1] = pres3;
        ylmd[s + 2] = 0.0;
        // l = 1, m = 0
        ylmd[s + 3] = 0.0;
        ylmd[s + 4] = 0.0;
        ylmd[s + 5] = pres3;
        // l = 1, m = 1
        ylmd[s + 6] = pres3;
        ylmd[s + 7] = 0.0;
        ylmd[s + 8] = 0.0;
        // clang-format on
    }
    if (this->i3sh.lmax < 2) return;

    // l = 2
    const Real pres15 = std::sqrt(15.0 / constants::PI4);
    const Real pres5h = std::sqrt(5.0 / constants::PI4) / 2.0;
    for (std::size_t i = 0; i < n; ++i)
    {
        std::size_t s = i * ldim + 4;
        const Real  x = rhat[3 * i];
        const Real  y = rhat[3 * i + 1];
        const Real  z = rhat[3 * i + 2];
        // clang-format off
        ylm[s    ] = pres15 * x * y;                 // l = 2, m = -2
        ylm[s + 1] = pres15 * y * z;                 // l = 2, m = -1
        ylm[s + 2] = pres5h * (3.0 * z * z - 1.0);   // l = 2, m = 0
        ylm[s + 3] = pres15 * x * z;                 // l = 2, m = 1
        ylm[s + 4] = pres15 * 0.5 * (x * x - y * y); // l = 2, m = 2
        s = 3 * i * ldim + 12;
        // l = 2, m = -2
        ylmd[s     ] = pres15 * y;
        ylmd[s +  1] = pres15 * x;
        ylmd[s +  2] = 0.0;
        // l = 2, m = -1
        ylmd[s +  3] = 0.0;
        ylmd[s +  4] = pres15 * z;
        ylmd[s +  5] = pres15 * y;
        // l = 2, m = 0
        ylmd[s +  6] = 0.0;
        ylmd[s +  7] = 0.0;
        ylmd[s +  8] = pres5h * 6.0 * z;
        // l = 2, m = 1
        ylmd[s +  9] = pres15 * z;
        ylmd[s + 10] = 0.0;
        ylmd[s + 11] = pres15 * x;
        // l = 2, m = 2
        ylmd[s + 12] = pres15 * x;
        ylmd[s + 13] = -pres15 * y;
        ylmd[s + 14] = 0.0;
        // clang-format on
    }
    if (this->i3sh.lmax < 3) return;

    // Zero all remaining entries for l > 2.
    for (std::size_t i = 0; i < n; ++i)
    {
        for (std::size_t ilm = 3 * 3; ilm < (std::size_t)ldim; ++ilm) ylm[i * ldim + ilm] = 0.0;
        for (std::size_t ilm = 3 * 3 * 3; ilm < 3 * (std::size_t)ldim; ++ilm)
        {
            ylmd[3 * i * ldim + ilm] = 0.0;
        }
    }

    // For l > 2 we use (a variant of) Clebsch-Gordan coefficients to recursively obtain the
    // remaining spherical harmonics.
    std::size_t l2 = 1;
    for (std::size_t l1 = 2; l1 < (std::size_t)i3sh.lmax; ++l1)
    {
        std::size_t       lmindx = i3sh.lookup[l1][l2];
        const std::size_t lsum = l1 + l2;
        for (std::size_t m1 = 0; m1 <= 2 * l1; ++m1)
        {
            for (std::size_t m2 = 0; m2 <= 2 * l2; ++m2)
            {
                for (Int ic = i3sh.indcg[lmindx]; ic < i3sh.indcg[lmindx + 1]; ++ic)
                {
                    const std::size_t lm = i3sh.js[ic];
                    if (lm + 1 > lsum * lsum && lm + 1 <= (lsum + 1) * (lsum + 1))
                    {
                        const Real c = i3sh.ylm3i[ic];
                        for (std::size_t i = 0; i < n; ++i)
                        {
                            const std::size_t s = i * ldim;
                            std::size_t       s1 = s + l1 * l1 + m1;
                            std::size_t       s2 = s + l2 * l2 + m2;
                            const std::size_t slm = 3 * (s + lm);
                            const Real        ylm1 = ylm[s1];
                            const Real        ylm2 = ylm[s2];
                            s1 *= 3;
                            s2 *= 3;
                            ylm[s + lm] += c * ylm1 * ylm2;
                            // clang-format off
                            ylmd[slm    ] += c * (ylmd[s1    ] * ylm2 + ylmd[s2    ] * ylm1);
                            ylmd[slm + 1] += c * (ylmd[s1 + 1] * ylm2 + ylmd[s2 + 1] * ylm1);
                            ylmd[slm + 2] += c * (ylmd[s1 + 2] * ylm2 + ylmd[s2 + 2] * ylm1);
                            // clang-format on
                        }
                    }
                }
                lmindx++;
            }
        }
    }

    return;
}

VASPML_NV_HOST_DEVICE
void math::SphericalHarmonics::computeSphericalHarmonicsAndGradients(const Vec1Real& rhat,
                                                                     const Vec1Real& rnorm,
                                                                     Vec1Real&       ylm,
                                                                     Vec1Real&       ylmd) const
{
#ifndef USE_NVIDIA_GPU
    VASPML_DEBUG_L1(
        if (rhat.size() != 3 * rnorm.size())
        {
            throw std::runtime_error("ERROR: Cannot compute spherical harmonics: input vectors "
                                     "do not match in size.");
        }
        if (ylm.size() != ldim * rnorm.size())
        {
            throw std::runtime_error("ERROR: Cannot compute spherical harmonics: incorrect size "
                                     "of output vector." + std::to_string(ylm.size()) +
                                     " " + std::to_string(ldim * rnorm.size()));
        }
        if (ylmd.size() != 3 * ldim * rnorm.size())
        {
            throw std::runtime_error("ERROR: Cannot compute spherical harmonics: incorrect size "
                                     "of derivative output vector.");
        }
    );
#endif

    // First, compute spherical harmonics and gradients for normalized vectors.
    computeSphericalHarmonicsAndGradientsNormalized(rhat, ylm, ylmd);

    // Finally, compute gradients with respect to unnormalized coordinates.
    for (std::size_t i = 0; i < rnorm.size(); ++i)
    {
        std::size_t const a = 3 * i;
        const Real        x = rhat[a];
        const Real        y = rhat[a + 1];
        const Real        z = rhat[a + 2];
        const Real        invrnorm = 1.0 / rnorm[i];
        std::size_t       s = a * ldim;
        for (std::size_t ilm = 0; ilm < (std::size_t)ldim; ++ilm)
        {
            const Real ylmdx = ylmd[s];
            const Real ylmdy = ylmd[s + 1];
            const Real ylmdz = ylmd[s + 2];
            const Real ylmdxa = ylmdx * x + ylmdy * y + ylmdz * z;
            // clang-format off
            ylmd[s    ] = (ylmdx - ylmdxa * x) * invrnorm;
            ylmd[s + 1] = (ylmdy - ylmdxa * y) * invrnorm;
            ylmd[s + 2] = (ylmdz - ylmdxa * z) * invrnorm;
            s += 3;
            // clang-format on
        }
    }

    return;
}

void math::SphericalHarmonics::integrate3SphericalHarmonics(Int                               lmax,
                                                            Integrate3SphericalHarmonicsData& data)
{
    VASPML_DEBUG_L1(
        if (data.ylm3.size() != 0 || data.ylm3i.size() != 0 || data.jl.size() != 0
            || data.js.size() != 0 || data.indcg.size() != 0 || data.lookup.size() != 0)
        {
            throw std::runtime_error("ERROR: Cannot integrate three spherical harmonics:"
                                     " data structure is not empty on entry to function.");
        }
    );

    using namespace constants;
    using namespace std;
    data.lmax = lmax;
    Vec1Real fn(4 * lmax + 2);
    math::factorial(fn);
    // Small helper function, if n is even return 1, else -1.
    auto fs = [](const Int n) { return n % 2 == 0 ? 1 : -1; };

    for (Int l1 = 0; l1 <= lmax; ++l1)
    {
        for (Int l2 = 0; l2 <= lmax; ++l2)
        {
            Real k2 = sqrt((2 * l1 + 1) * (2 * l2 + 1));
            for (Int m1 = -l1; m1 <= l1; ++m1)
            {
                for (Int m2 = -l2; m2 <= l2; ++m2)
                {
                    data.indcg.push_back(data.ylm3.size());

                    const Int n1 = abs(m1);
                    const Int n2 = abs(m2);
                    const Int s1 = m1 < 0 ? 1 : 0;
                    const Int s2 = m2 < 0 ? 1 : 0;
                    const Int t1 = m1 == 0 ? 1 : 0;
                    const Int t2 = m2 == 0 ? 1 : 0;

                    // There are potentially two non-zero m3 values.
                    Int m3;
                    Int m3p;
                    Int nm3;
                    if (m1 * m2 < 0)
                    {
                        m3 = -n1 - n2;
                        m3p = -abs(n1 - n2);
                        if (m3p == 0) nm3 = 1;
                        else nm3 = 2;
                    }
                    else if (m1 * m2 == 0)
                    {
                        m3 = m1 + m2;
                        m3p = 0;
                        nm3 = 1;
                    }
                    else
                    {
                        m3 = n1 + n2;
                        m3p = abs(n1 - n2);
                        nm3 = 2;
                    }

                    do {
                        const Int n3 = abs(m3);
                        const Int s3 = m3 < 0 ? 1 : 0;
                        const Int t3 = m3 == 0 ? 1 : 0;

                        const Real q1 = 0.5 * k2 * fs(n3 + (s1 + s2 + s3) / 2);
                        const Real q2 = 1.0 / pow(sqrt(2.0), 1 + t1 + t2 + t3);

                        for (Int l3 = abs(l1 - l2); l3 <= l1 + l2; l3 += 2)
                        {
                            if (n3 > l3) continue;
                            Real t = 0.0;
                            // clang-format off
                            if (n1 + n2 == -n3) t += clebschGordanM0(l1, l2, l3, fn);
                            if (n1 + n2 ==  n3) t += clebschGordan(l1, l2, l3, n1, n2, n3, fn)
                                                   * fs(n3 + s3);
                            if (n1 - n2 == -n3) t += clebschGordan(l1, l2, l3, n1, -n2, -n3, fn)
                                                   * fs(n2 + s2);
                            if (n1 - n2 ==  n3) t += clebschGordan(l1, l2, l3, -n1, n2, -n3, fn)
                                                   * fs(n1 + s1);
                            // clang-format on
                            const Real t0 = clebschGordanM0(l1, l2, l3, fn);
                            const Real tmp = SQRT_PI * sqrt(2 * l3 + 1);
                            data.ylm3.push_back(q1 * q2 * t * t0 / tmp);
                            if (abs(t0) < 1.0E-15) data.ylm3i.push_back(0.0);
                            else data.ylm3i.push_back(t * q2 * tmp / (q1 * t0));
                            data.jl.push_back(l3);
                            data.js.push_back(l3 * (l3 + 1) + m3);
                        }
                        nm3--;
                        m3 = m3p;
                    }
                    while (nm3 > 0);
                }
            }
        }
    }

    data.indcg.push_back(data.ylm3.size());
    Int lmind = 0;
    for (Int l1 = 0; l1 <= lmax; ++l1)
    {
        data.lookup.push_back(Vec1Int());
        for (Int l2 = 0; l2 <= lmax; ++l2)
        {
            data.lookup.back().push_back(lmind);
            lmind += (2 * l1 + 1) * (2 * l2 + 1);
        }
    }

    return;
}

std::size_t math::SphericalHarmonics::get_lmax(void) const
{
    return i3sh.lmax;
}

std::size_t math::SphericalHarmonics::get_ldim(void) const
{
    return ldim;
}

const math::SphericalHarmonics::Integrate3SphericalHarmonicsData&
math::SphericalHarmonics::get_i3sh() const
{
    return i3sh;
}
