#ifndef __NEC__
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MODULE MultiType

#include "types.hpp"

#include "boost_helpers.hpp"

#include <boost/test/unit_test.hpp>

using namespace vaspml;

BOOST_AUTO_TEST_SUITE(UnitTests)

BOOST_AUTO_TEST_CASE(NonSharedPtrNonConst_CorrectValuesAndAddresses)
{
    MultiType a = 1;

    Int&       aref = a.get<Int>();
    const Int& acref = a.cget<Int>();
    a = 2;

    // Check value.
    REQUIRE_EQUAL_MESSAGE(aref, 2, "value aref");
    REQUIRE_EQUAL_MESSAGE(acref, 2, "value acref");
    // Check address.
    REQUIRE_EQUAL_MESSAGE(&aref, &(std::get<Int>(a)), "address aref");
    REQUIRE_EQUAL_MESSAGE(&acref, &(std::get<Int>(a)), "address acref");
}

BOOST_AUTO_TEST_CASE(NonSharedPtrConst_CorrectValuesAndAddresses)
{
    const MultiType a = 1;

    const Int& acref = a.cget<Int>();

    // Check value.
    REQUIRE_EQUAL_MESSAGE(acref, 1, "value acref");
    // Check address.
    REQUIRE_EQUAL_MESSAGE(&acref, &(std::get<Int>(a)), "address acref");
}

BOOST_AUTO_TEST_CASE(SharedPtrNonConst_CorrectValuesAndAddresses)
{
    MultiType v = std::make_shared<Vec1Int>(Vec1Int({1, 2, 3}));

    ShVec1Int         v1 = v.get<ShVec1Int>();
    ShVec1Int&        v2 = v.get<ShVec1Int>();
    Vec1Int&          v3 = v.dget<ShVec1Int>();
    ShCVec1Int        c1 = v.cget<ShVec1Int>();
    const ShCVec1Int& c2 = v.cget<ShVec1Int>();
    const Vec1Int&    c3 = v.dcget<ShVec1Int>();

    // Check values.
    REQUIRE_EQUAL_COLLECTIONS(*std::get<ShVec1Int>(v), *v1, "values v1");
    REQUIRE_EQUAL_COLLECTIONS(*std::get<ShVec1Int>(v), *v2, "values v2");
    REQUIRE_EQUAL_COLLECTIONS(*std::get<ShVec1Int>(v), v3, "values v3");
    REQUIRE_EQUAL_COLLECTIONS(*std::get<ShVec1Int>(v), *c1, "values c1");
    REQUIRE_EQUAL_COLLECTIONS(*std::get<ShVec1Int>(v), *c2, "values c2");
    REQUIRE_EQUAL_COLLECTIONS(*std::get<ShVec1Int>(v), c3, "values c3");
    // Check vector address.
    REQUIRE_EQUAL_MESSAGE(std::get<ShVec1Int>(v), v1, "vector address v1");
    REQUIRE_EQUAL_MESSAGE(std::get<ShVec1Int>(v), v2, "vector address v2");
    REQUIRE_EQUAL_MESSAGE(std::get<ShVec1Int>(v).get(), &v3, "vector address v3");
    REQUIRE_EQUAL_MESSAGE(std::get<ShVec1Int>(v), c1, "vector address c1");
    REQUIRE_EQUAL_MESSAGE(std::get<ShVec1Int>(v), c2, "vector address c2");
    REQUIRE_EQUAL_MESSAGE(std::get<ShVec1Int>(v).get(), &c3, "vector address c3");
    // Check data starting address.
    REQUIRE_EQUAL_MESSAGE(std::get<ShVec1Int>(v)->data(), v1->data(), "data staring address v1");
    REQUIRE_EQUAL_MESSAGE(std::get<ShVec1Int>(v)->data(), v2->data(), "data staring address v2");
    REQUIRE_EQUAL_MESSAGE(std::get<ShVec1Int>(v)->data(), v3.data(), "data staring address v3");
    REQUIRE_EQUAL_MESSAGE(std::get<ShVec1Int>(v)->data(), c1->data(), "data staring address c1");
    REQUIRE_EQUAL_MESSAGE(std::get<ShVec1Int>(v)->data(), c2->data(), "data staring address c2");
    REQUIRE_EQUAL_MESSAGE(std::get<ShVec1Int>(v)->data(), c3.data(), "data staring address c3");
}

BOOST_AUTO_TEST_CASE(SharedPtrConst_CorrectValuesAndAddresses)
{
    const MultiType v = std::make_shared<Vec1Int>(Vec1Int({1, 2, 3}));

    ShCVec1Int        c1 = v.cget<ShVec1Int>();
    const ShCVec1Int& c2 = v.cget<ShVec1Int>();
    const Vec1Int&    c3 = v.dcget<ShVec1Int>();

    // Check values.
    REQUIRE_EQUAL_COLLECTIONS(*std::get<ShVec1Int>(v), *c1, "values c1");
    REQUIRE_EQUAL_COLLECTIONS(*std::get<ShVec1Int>(v), *c2, "values c2");
    REQUIRE_EQUAL_COLLECTIONS(*std::get<ShVec1Int>(v), c3, "values c3");
    // Check vector address.
    REQUIRE_EQUAL_MESSAGE(std::get<ShVec1Int>(v), c1, "vector address c1");
    REQUIRE_EQUAL_MESSAGE(std::get<ShVec1Int>(v), c2, "vector address c2");
    REQUIRE_EQUAL_MESSAGE(std::get<ShVec1Int>(v).get(), &c3, "vector address c3");
    // Check data starting address.
    REQUIRE_EQUAL_MESSAGE(std::get<ShVec1Int>(v)->data(), c1->data(), "data starting address c1");
    REQUIRE_EQUAL_MESSAGE(std::get<ShVec1Int>(v)->data(), c2->data(), "data starting address c2");
    REQUIRE_EQUAL_MESSAGE(std::get<ShVec1Int>(v)->data(), c3.data(), "data starting address c3");
}

BOOST_AUTO_TEST_CASE(AssignExisting_ObjectAliased)
{
    ShVec1Int a = std::make_shared<Vec1Int>(Vec1Int{1, 2, 3, 4, 5});
    REQUIRE_EQUAL_MESSAGE(a.use_count(), 1, "a.use_count() before");

    ShVec1Int b = assignOrMakeShared(a);
    REQUIRE_EQUAL_MESSAGE(a.use_count(), 2, "a.use_count() after");
    REQUIRE_EQUAL_MESSAGE(b.use_count(), 2, "b.use_count() after");
    REQUIRE_EQUAL_MESSAGE(a, b, "vector address");
    REQUIRE_EQUAL_MESSAGE(a->data(), b->data(), "data starting address");
    Vec1Int& va = *a;
    Vec1Int& vb = *b;
    REQUIRE_EQUAL_COLLECTIONS(va, vb, "values");
    va[1] = 99;
    va[4] = 199;
    REQUIRE_EQUAL_COLLECTIONS(va, vb, "values modified");
}

BOOST_AUTO_TEST_CASE(AssignExistingReference_ObjectAliased)
{
    ShVec1Int  a = std::make_shared<Vec1Int>(Vec1Int{1, 2, 3, 4, 5});
    ShVec1Int& c = a;
    REQUIRE_EQUAL_MESSAGE(a.use_count(), 1, "a.use_count() before");

    ShVec1Int b = assignOrMakeShared(c);
    REQUIRE_EQUAL_MESSAGE(a.use_count(), 2, "a.use_count() after");
    REQUIRE_EQUAL_MESSAGE(b.use_count(), 2, "b.use_count() after");
    REQUIRE_EQUAL_MESSAGE(a, b, "vector address");
    REQUIRE_EQUAL_MESSAGE(a->data(), b->data(), "data starting address");
    Vec1Int& va = *a;
    Vec1Int& vb = *b;
    REQUIRE_EQUAL_COLLECTIONS(va, vb, "values");
    va[1] = 99;
    va[4] = 199;
    REQUIRE_EQUAL_COLLECTIONS(va, vb, "values modified");
}

BOOST_AUTO_TEST_CASE(AssignExistingConstReference_ObjectAliased)
{
    ShVec1Int        a = std::make_shared<Vec1Int>(Vec1Int{1, 2, 3, 4, 5});
    const ShVec1Int& c = a;
    REQUIRE_EQUAL_MESSAGE(a.use_count(), 1, "a.use_count() before");

    ShVec1Int b = assignOrMakeShared(c);
    REQUIRE_EQUAL_MESSAGE(a.use_count(), 2, "a.use_count() after");
    REQUIRE_EQUAL_MESSAGE(b.use_count(), 2, "b.use_count() after");
    REQUIRE_EQUAL_MESSAGE(a, b, "vector address");
    REQUIRE_EQUAL_MESSAGE(a->data(), b->data(), "data starting address");
    Vec1Int& va = *a;
    Vec1Int& vb = *b;
    REQUIRE_EQUAL_COLLECTIONS(va, vb, "values");
    va[1] = 99;
    va[4] = 199;
    REQUIRE_EQUAL_COLLECTIONS(va, vb, "values modified");
}

BOOST_AUTO_TEST_CASE(AssignNullptr_NewObjectCreated)
{
    ShVec1Int a = std::make_shared<Vec1Int>(Vec1Int{1, 2, 3, 4, 5});
    REQUIRE_EQUAL_MESSAGE(a.use_count(), 1, "a.use_count() before");
    a = nullptr;
    REQUIRE_EQUAL_MESSAGE(a.use_count(), 0, "a.use_count() assigned nullptr");

    ShVec1Int b = assignOrMakeShared(a);
    REQUIRE_EQUAL_MESSAGE(a.use_count(), 0, "a.use_count() after");
    REQUIRE_EQUAL_MESSAGE(b.use_count(), 1, "b.use_count() after");
}

BOOST_AUTO_TEST_SUITE_END()
