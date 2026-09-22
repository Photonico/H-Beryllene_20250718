#ifndef __NEC__
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MODULE Lattice

#include "Lattice.hpp"
#include "boost_helpers.hpp"

#include <boost/test/unit_test.hpp>

using namespace vaspml;

BOOST_AUTO_TEST_SUITE(UnitTests)

BOOST_AUTO_TEST_CASE(InverseOfInverseLattice_Identity)
{
    Lattice lattice(1.0, 0.3, 0.2, 0.4, 1.5, 0.1, 0.1, 0.05, 2.0);
    Lattice inverse = lattice.computeInverse();

    for (std::size_t i = 0; i < 3; ++i)
    {
        for (std::size_t j = 0; j < 3; ++j)
        {
            Real mij = 0.0;
            for (std::size_t k = 0; k < 3; ++k)
            {
                mij += lattice.component(i, k) * inverse.component(k, j);
            }
            if (i == j) BOOST_REQUIRE_SMALL(mij - 1.0, 1.0E-15);
            else BOOST_REQUIRE_SMALL(mij, 1.0E-15);
        }
    }
}

BOOST_AUTO_TEST_CASE(DirectCartesianDirect_OriginalVectorRestored)
{
    Lattice lattice;
    lattice.set_components(Vec1Real{1.0, 0.3, 0.2, 0.4, 1.5, 0.1, 0.1, 0.05, 2.0});
    Lattice        inverse = lattice.computeInverse();
    const Vec1Real dvec{0.3, 0.4, 0.8};

    Vec1Real result = lattice.timesVector(dvec);
    result = inverse.timesVector(result);

    REQUIRE_CLOSE_COLLECTIONS(dvec, result, 1.0E-15, "dvec");
}

BOOST_AUTO_TEST_SUITE_END()
