#include "EquationSolver.cpp"

double f(double x)
{
    double L = 10;
    double r = 1;
    double V = 12.4;
    return L * (0.5*M_PI*std::pow(r,2) - std::pow(r,2)*std::asin(x/r) - x*std::sqrt(r*r-x*x)) - V;
}

double df(double x)
{
    double L = 10;
    double r = 1;
    double V = 12.4;
    return L * (-2*x*std::asin(x/r) - x*x/r/std::sqrt(1 - std::pow(x/r,2)) - std::sqrt(r*r-x*x)
        + x*x/std::sqrt(r*r-x*x));
}

int main(int argc , char* argv[])
{
    class bisection solver_bi(f, 0, 1);
    class Newton solver_Newton(f,df,0.5);
    class secant solver_secant(f,0.25,0.75);
    double x_bi = solver_bi.solve();    double fx_bi = f(x_bi);
    double x_Newton = solver_Newton.solve();    double fx_Newton = f(x_Newton);
    double x_secant = solver_secant.solve();    double fx_secant = f(x_secant);
    std::cout << "bisection: " << x_bi << "; fx: " << fx_bi << std::endl;
    std::cout << "Newton: " << x_Newton << "; fx: " << fx_Newton << std::endl;
    std::cout << "secant: " << x_secant << "; fx: " << fx_secant << std::endl;
    return 0;
}