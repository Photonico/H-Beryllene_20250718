#ifndef __NEC__
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MODULE LinearAlgebra

#include "Linalg.hpp"
#include "TestCase_Linalg.hpp"
#include "boost_helpers.hpp"

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

using namespace vaspml;
namespace bdata = boost::unit_test::data;

TestCaseContainer<TestCaseLinalg_l2Norm>                container_l2Norm;
TestCaseContainer<TestCaseLinalg_scaleVector>           container_scaleVector;
TestCaseContainer<TestCaseLinalg_matMul>                container_matMul;
TestCaseContainer<TestCaseLinalg_scaleVectorPlusVector> container_scaleVectorPlusVector;
TestCaseContainer<TestCaseLinalg_dotProduct>            container_dotProduct;
linalg::LinalgContext                                   context;

BOOST_AUTO_TEST_SUITE(UnitTests)
BOOST_DATA_TEST_CASE(l2Norm_Vector, bdata::make(container_l2Norm.testCases), testCase)
{

    Vec1Real norm;
    norm.push_back(linalg::l2Norm(testCase.vector[0], (Int)testCase.vector[0].size(), context));
    norm.push_back(linalg::l2Norm(testCase.vector[1], (Int)testCase.vector[1].size(), context));
    //norm.push_back( linalg::l2Norm( testCase.vector[1] ) );
    REQUIRE_CLOSE_COLLECTIONS(norm, testCase.norm, testCase.tolerance, "linalg::l2Norm");
}

BOOST_DATA_TEST_CASE(scaleVector, bdata::make(container_scaleVector.testCases), testCase)
{

    Vec1Real vector = testCase.vector;
    linalg::scaleVector(vector, testCase.scalar, vector.size(), context);
    REQUIRE_CLOSE_COLLECTIONS(vector, testCase.target, testCase.tolerance, "linalg::scaleVector");
}

BOOST_DATA_TEST_CASE(matMul_rectangularMatrix, bdata::make(container_matMul.testCases), testCase)
{

    Vec1Real MatrixC;
    MatrixC.resize(20);
    linalg::matMul(testCase.n,
                   testCase.k,
                   testCase.m,
                   testCase.MatrixA,
                   testCase.MatrixB,
                   MatrixC,
                   context);
    REQUIRE_CLOSE_COLLECTIONS(MatrixC, testCase.MatrixC, testCase.tolerance, "linalg::matMul");
}
BOOST_AUTO_TEST_SUITE_END()

BOOST_DATA_TEST_CASE(scaleVectorPlusVector,
                     bdata::make(container_scaleVectorPlusVector.testCases),
                     testCase)
{

    Vec1Real testVector = testCase.addVector;
    linalg::scaleVectorPlusVector(testCase.scalar,
                                  testCase.scaleVector,
                                  testVector,
                                  testVector.size(),
                                  context);
    REQUIRE_CLOSE_COLLECTIONS(testVector,
                              testCase.target,
                              testCase.tolerance,
                              "linalg::scaleVectorPlusVector");
}

BOOST_DATA_TEST_CASE(dotProduct, bdata::make(container_dotProduct.testCases), testCase)
{

    Real norm =
        linalg::dotProduct(testCase.vectorA, testCase.vectorB, testCase.vectorB.size(), context);
    BOOST_REQUIRE_SMALL(norm - testCase.target, testCase.tolerance);
}
