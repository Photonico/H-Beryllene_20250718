#ifndef TESTCASE_HPP
#define TESTCASE_HPP

#include <ostream>
#include <string>

namespace vaspml
{

struct TestCase
{
    std::string name;

    TestCase(std::string name);
};

std::ostream& operator<<(std::ostream& os, vaspml::TestCase const& testCase);

} //namespace vaspml

#endif
