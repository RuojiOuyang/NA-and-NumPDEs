#include "EquationSolver.cpp"

double fA1(double x)
{ 
    return 1 / x - x * x;
}

int main(int argc, char** argv)
{
    class bisection A(fA1, 1, 2);
    A.test();
}