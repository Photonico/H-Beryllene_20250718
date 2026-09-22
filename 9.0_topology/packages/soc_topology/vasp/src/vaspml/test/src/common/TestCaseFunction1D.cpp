#include "TestCaseFunction1D.hpp"

using namespace vaspml;

TestCaseFunction1D::TestCaseFunction1D(std::string name) :
    TestCase(name),
    tolerance(defaultTolerance),
    dtolerance(defaultTolerance)
{}
