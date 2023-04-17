#ifndef __TEST_H__
#define __TEST_H__

#include <cppunit/extensions/HelperMacros.h>
#include "PolyInterp.h"

class PolyTest : public CppUnit::TestFixture 
{
    CPPUNIT_TEST_SUITE( PolyTest );

    CPPUNIT_TEST( polyEquiv1 ); ///< test p1 == p2
    CPPUNIT_TEST( polyEquiv2 ); ///< test p1 == p1

    CPPUNIT_TEST( polyAssign1 ); ///< test p2 = p1
    CPPUNIT_TEST( polyAssign2 ); ///< test p1 = c(real number)

    CPPUNIT_TEST( polyAdd1 ); ///< test p1 + p2
    CPPUNIT_TEST( polyAdd2 ); ///< test p1 + c
    CPPUNIT_TEST( polyAdd3 ); ///< test c + p1

    CPPUNIT_TEST( polyMinus1 ); ///< test -p1
    CPPUNIT_TEST( polyMinus2 ); ///< test p1 - p2
    CPPUNIT_TEST( polyMinus3 ); ///< test p1 - c
    CPPUNIT_TEST( polyMinus4 ); ///< test c - p1

    CPPUNIT_TEST( polyMulti1 ); ///< test p1 * p2
    CPPUNIT_TEST( polyMulti2 ); ///< test p1 * c
    CPPUNIT_TEST( polyMulti3 ); ///< test c * p1

    CPPUNIT_TEST( polyDiff1 ); ///< test p1'
    CPPUNIT_TEST( polyDiff2 ); ///< test p1''

    CPPUNIT_TEST( polyDegree ); ///< test \partial(p1)

    CPPUNIT_TEST( polyEvaluation ); ///< test p1(x_0)

    CPPUNIT_TEST_SUITE_END();

private:
    Polynomial p1;
    Polynomial p2;

public:
    void setUp();  

    void polyEquiv1();
    void polyEquiv2();

    void polyAssign1();
    void polyAssign2();

    void polyAdd1();
    void polyAdd2();
    void polyAdd3();

    void polyMinus1();
    void polyMinus2();
    void polyMinus3();
    void polyMinus4();

    void polyMulti1();
    void polyMulti2();
    void polyMulti3();

    void polyDiff1();
    void polyDiff2();

    void polyDegree();

    void polyEvaluation();
};

class InterpTest : public CppUnit::TestFixture 
{
    CPPUNIT_TEST_SUITE(InterpTest);
    CPPUNIT_TEST( Newton_overwrite );
    CPPUNIT_TEST( Newton_addition );
    CPPUNIT_TEST( Neville_Aitken );
    CPPUNIT_TEST_SUITE_END();

private:
    InterpConditions cond1;
    InterpConditions cond2;
    InterpConditions cond3;
    InterpConditions cond4;

public:
    InterpTest();
    void Newton_overwrite();
    void Newton_addition();
    void Neville_Aitken();
};
#endif