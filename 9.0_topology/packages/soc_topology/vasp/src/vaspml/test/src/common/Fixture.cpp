#include "Fixture.hpp"

#include <boost/test/unit_test.hpp>

using namespace vaspml;

Fixture::Fixture(std::string contentDescription)
{
    BOOST_TEST_MESSAGE("Loading fixture \"" + contentDescription + "\".");
}
