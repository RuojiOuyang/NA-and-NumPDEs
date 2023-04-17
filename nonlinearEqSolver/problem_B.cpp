#include "EquationSolver.cpp"

double f1(double x)
{ 
    return 1 / x - std::tan(x);
}

double f2(double x)
{
    return 1 / x - std::pow(2 , x);
}

double f3(double x)
{
    return std::pow(2 , -x) + std::exp(x) + 2 * std::cos(x) - 6;
}

double f4(double x)
{
    return (std::pow(x , 3) + 4 * std::pow(x , 2) + 3 * x + 5) / (2 * std::pow(x , 3) - 9 * std::pow(x , 2) + 18 * x - 2);
}

double f4_denominator(double x)
{
    return (2 * std::pow(x , 3) - 9 * std::pow(x , 2) + 18 * x - 2);
}

int main(int argc, char** argv)
{
    class bisection solver1(f1, 0, M_PI/2);
    class bisection solver2(f2 , 0 , 1);
    class bisection solver3(f3 , 1 , 3);
    class bisection solver4(f4 , 0 , 4);
    double x1 = solver1.solve(); double fx1 = f1(x1);
    double x2 = solver2.solve(); double fx2 = f2(x2);
    double x3 = solver3.solve(); double fx3 = f3(x3);
    double x4 = solver4.solve(); double fx4 = f4(x4);
    std::cout << x1 << " " << fx1 << std::endl;
    std::cout << x2 << " " << fx2 << std::endl;
    std::cout << x3 << " " << fx3 << std::endl;
    std::cout << x4 << " " << fx4 << std::endl;

    class bisection solver5(f4_denominator, 0 , 4);
    double x5 = solver5.solve(); double fx5 = f4_denominator(x4);
    std::cout << x5 << " " << fx5 << std::endl;
    return 0;
}