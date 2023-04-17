#include <cppunit/extensions/HelperMacros.h>
#include <iostream>
#include <cstring>

#include "test.h"

InterpTest::InterpTest()
{
    double x[3] = {0 , 1 , 3};
    double_array f1 = {1};
    double_array f2 = {2 , -1};
    double_array f3 = {0 , 0};

    ///< cond1
    (cond1.x).assign(x , x+3);
    (cond1.fx).push_back(f1);
    (cond1.fx).push_back(f2);
    (cond1.fx).push_back(f3);

    ///< cond2
    (cond2.x).assign(x , x+2);
    (cond2.fx).push_back(f1);
    (cond2.fx).push_back(f2);

    ///< cond3
    (cond3.x).assign(x+2, x+3);
    (cond3.fx).push_back(f3);

    ///< cond4
    double_array x_ = {0, 1, 2, 3};
    double2_array fx_ = {{6},{-3},{-6},{9}};
    cond4.x = x_;
    cond4.fx = fx_;
}

void InterpTest::Newton_overwrite()
{
    NewtonInterp myInterp;
    double_array coeff = {1, 49/12.0, -155/36.0, 49/36.0, -5/36.0};
    Polynomial expected(coeff);
    myInterp.Newton_overwrite(cond1);
    CPPUNIT_ASSERT(myInterp.get_Polynomial() == expected);
}

void InterpTest::Newton_addition()
{
    NewtonInterp myInterp;
    double_array coeff = {1, 49/12.0, -155/36.0, 49/36.0, -5/36.0};
    Polynomial expected(coeff);
    myInterp.Newton_overwrite(cond2);
    myInterp.Newton_addition(cond3);
    CPPUNIT_ASSERT(myInterp.get_Polynomial() == expected);
}

void InterpTest::Neville_Aitken()
{
    CPPUNIT_ASSERT(NewtonInterp::Neville_Aitken(cond4, 1.5) == -6);
}