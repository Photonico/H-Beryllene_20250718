#ifndef __NEC__
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MODULE utils

#include <boost/test/unit_test.hpp>

#include "boost_helpers.hpp"
#include "utils.hpp"

using namespace vaspml;

BOOST_AUTO_TEST_SUITE(UnitTests)

BOOST_AUTO_TEST_CASE(TrimDifferentStrings_TrimmedStrings)
{
    BOOST_REQUIRE_EQUAL(string_tools::ltrim("abc def ghi"), "abc def ghi");
    BOOST_REQUIRE_EQUAL(string_tools::ltrim("  abc def ghi"), "abc def ghi");
    BOOST_REQUIRE_EQUAL(string_tools::ltrim("abc def ghi  "), "abc def ghi  ");
    BOOST_REQUIRE_EQUAL(string_tools::ltrim("  abc def ghi  "), "abc def ghi  ");

    BOOST_REQUIRE_EQUAL(string_tools::rtrim("abc def ghi"), "abc def ghi");
    BOOST_REQUIRE_EQUAL(string_tools::rtrim("  abc def ghi"), "  abc def ghi");
    BOOST_REQUIRE_EQUAL(string_tools::rtrim("abc def ghi  "), "abc def ghi");
    BOOST_REQUIRE_EQUAL(string_tools::rtrim("  abc def ghi  "), "  abc def ghi");

    BOOST_REQUIRE_EQUAL(string_tools::trim("abc def ghi"), "abc def ghi");
    BOOST_REQUIRE_EQUAL(string_tools::trim("  abc def ghi"), "abc def ghi");
    BOOST_REQUIRE_EQUAL(string_tools::trim("abc def ghi  "), "abc def ghi");
    BOOST_REQUIRE_EQUAL(string_tools::trim("  abc def ghi  "), "abc def ghi");
}

BOOST_AUTO_TEST_CASE(CreateStringsFromPrintfStyle_CorrectSizesAndContents)
{
    BOOST_REQUIRE_EQUAL(str("%04d", 44), "0044");
    BOOST_REQUIRE_EQUAL(str("%04d", 44).size(), 4);
    BOOST_REQUIRE_EQUAL(str("%4s", "ab"), "  ab");
    BOOST_REQUIRE_EQUAL(str("%4s", "ab").size(), 4);
    BOOST_REQUIRE_EQUAL(str("%-4s", "ab"), "ab  ");
    BOOST_REQUIRE_EQUAL(str(("%04d" + std::string(STR_MAX - 5, '?')).c_str(), 44),
                        "0044" + std::string(STR_MAX - 5, '?'));
}

BOOST_AUTO_TEST_CASE(SplitString_CheckCorrectParts)
{
    Vec1String result = string_tools::splitString("bacadae", "a");
    Vec1String expect = {"b", "c", "d", "e"};
    REQUIRE_EQUAL_COLLECTIONS(result, expect, "bacadae");

    result = string_tools::splitString("aabaacaadaaeaa", "aa", false);
    expect = {"", "b", "c", "d", "e", ""};
    REQUIRE_EQUAL_COLLECTIONS(result, expect, "aabaacaadaaeaa");

    result = string_tools::splitString("aabaacaadaaeaa", "aa");
    expect = {"b", "c", "d", "e"};
    REQUIRE_EQUAL_COLLECTIONS(result, expect, "aabaacaadaaeaa");

    result = string_tools::splitString("a ; b ; c ; d ; e", ";");
    expect = {"a ", " b ", " c ", " d ", " e"};
    REQUIRE_EQUAL_COLLECTIONS(result, expect, "a ; b ; c ; d ; e");

    result = string_tools::splitString("abacadaea", "a");
    expect = {"b", "c", "d", "e"};
    REQUIRE_EQUAL_COLLECTIONS(result, expect, "abacadaea");
}

BOOST_AUTO_TEST_SUITE_END()
