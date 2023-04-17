#include "EquationSolver.cpp"

double f(double x)
{
    return tan(x) - x;
}

double df(double x)
{
    return pow(cos(x) , -2) - 1;
}

int main(int argc, char** argv)
{
    class Newton solver1(f , df , 4.5);
    class Newton solver2(f , df , 7.7);
    double x1 = solver1.solve(); double fx1 = f(x1);
    double x2 = solver2.solve(); double fx2 = f(x2);
    std::cout << "root near 4.5: " << x1 << "; fx:" << fx1 << std::endl;
    std::cout << "root near 7.7: " << x2 << "; fx:" << fx2 << std::endl;
    return 0;
}