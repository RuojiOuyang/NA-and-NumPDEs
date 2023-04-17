#include <cppunit/extensions/HelperMacros.h>
#include <iostream>
#include <cstring>

#include "test.h"

void PolyTest::setUp()
{
    double_array coeff1 = {3 , 1 , 1};
    double_array coeff2 = {4 , 3 , 1};
    p1.set_coefficient(coeff1);
    p2.set_coefficient(coeff2);
}

void PolyTest::polyEquiv1()
{
    CPPUNIT_ASSERT((p1 == p2) == false);
}

void PolyTest::polyEquiv2()
{
    CPPUNIT_ASSERT((p1 == p1) == true);
}

void PolyTest::polyAssign1()
{
    p2 = p1;
    CPPUNIT_ASSERT((p2 == p1) == true);
}

void PolyTest::polyAssign2()
{
    p1 = 1;
    double_array coeff = {1};
    Polynomial expected(coeff);
    CPPUNIT_ASSERT(p1 == expected);
}

void PolyTest::polyAdd1()
{
    double_array coeff = {7 , 4 , 2};
    Polynomial expected(coeff);
    CPPUNIT_ASSERT((p1 + p2) == expected);
}

void PolyTest::polyAdd2()
{
    double_array coeff = {7 , 1 , 1};
    Polynomial expected(coeff);
    CPPUNIT_ASSERT((p1 + (double)4) == expected);
}

void PolyTest::polyAdd3()
{
    double_array coeff = {7 , 1 , 1};
    Polynomial expected(coeff);
    CPPUNIT_ASSERT(((double)4 + p1) == expected);
}

void PolyTest::polyMinus1()
{
    double_array coeff = {-3 , -1 , -1};
    Polynomial expected(coeff);
    CPPUNIT_ASSERT((-p1) == expected);
}

void PolyTest::polyMinus2()
{
    double_array coeff = {-1 , -2};
    Polynomial expected(coeff);
    CPPUNIT_ASSERT((p1 - p2) == expected);
}

void PolyTest::polyMinus3()
{
    double_array coeff = {0 , 1 , 1};
    Polynomial expected(coeff);
    CPPUNIT_ASSERT((p1 - (double)3) == expected);
}

void PolyTest::polyMinus4()
{
    double_array coeff = {0 , -1 , -1};
    Polynomial expected(coeff);
    CPPUNIT_ASSERT(((double)3 - p1) == expected);
}

void PolyTest::polyMulti1()
{
    double_array coeff = {12, 13, 10, 4, 1};
    Polynomial expected(coeff);
    CPPUNIT_ASSERT((p1 * p2) == expected);
}

void PolyTest::polyMulti2()
{
    double_array coeff = {9 , 3 , 3};
    Polynomial expected(coeff);
    CPPUNIT_ASSERT((p1 * (double)3) == expected);
}

void PolyTest::polyMulti3()
{
    double_array coeff = {9 , 3 , 3};
    Polynomial expected(coeff);
    CPPUNIT_ASSERT(((double)3 * p1) == expected);
}

void PolyTest::polyDiff1()
{
    double_array coeff = {1 , 2};
    Polynomial expected(coeff);
    CPPUNIT_ASSERT(p1.diff() == expected);
}

void PolyTest::polyDiff2()
{
    double_array coeff = {2};
    Polynomial expected(coeff);
    CPPUNIT_ASSERT(p1.diff(2) == expected);
}

void PolyTest::polyDegree()
{
    CPPUNIT_ASSERT(p1.get_degree() == 2);
}

void PolyTest::polyEvaluation()
{
    CPPUNIT_ASSERT(p1(3) == 15);
}