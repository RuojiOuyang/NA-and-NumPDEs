#include "EquationSolver.cpp"

double l = 89;
double h = 49;
double D = 55;
double beta = 11.5 * M_PI / 180;   

double f(double x)
{
    double A = l * std::sin(beta);
    double B = l * std::cos(beta);
    double C = (h + 0.5 * D) * std::sin(beta) - 0.5 * D * std::tan(beta);
    double E = (h + 0.5 * D) * std::cos(beta) - 0.5 * D;
    return (A*std::sin(x)*std::cos(x)+B*std::pow(std::sin(x),2)-C*std::cos(x)-E*std::sin(x));
}

double df(double x)
{
    double A = l * std::sin(beta);
    double B = l * std::cos(beta);
    double C = (h + 0.5 * D) * std::sin(beta) - 0.5 * D * std::tan(beta);
    double E = (h + 0.5 * D) * std::cos(beta) - 0.5 * D;
    return (A*std::cos(2*x)+B*std::sin(2*x)+C*std::sin(x)-E*std::cos(x));
}

int main(int argc , char* argv[])
{
    ///--------------(a)---------------
    class Newton solver(f, df, 33*M_PI/180);
    double alpha = solver.solve(); double falpha = f(alpha);
    std::cout << "(a)"<<std::endl;
    std::cout << "alpha: " << alpha*180/M_PI << "; f(alpha): " << falpha << std::endl;

    ///--------------(b)---------------
    D = 30;
    alpha = solver.solve(); falpha = f(alpha);
    std::cout << "(b)"<<std::endl;
    std::cout << "alpha: " << alpha*180/M_PI << "; f(alpha): " << falpha << std::endl;

    ///--------------(c)---------------
    double initial_points[11][2] = {{100,200},{200,300},{300,400},{400,500},{500,600},{600,700},{700,800},{800,900},{900,1000},{1000,2000},{10000,20000}};
    class secant solver_secant(f, 0, 0);
    std::cout << "(c)"<<std::endl;
    for(int i = 0 ; i < 11 ; i ++)
    {
        double alpha0 = initial_points[i][0]* M_PI / 180; 
        double alpha1 = initial_points[i][1]* M_PI / 180; 
        solver_secant.set_initial(alpha0, alpha1);
        alpha = solver_secant.solve(); falpha = f(alpha);
        printf("%d. initial points: %.0f, %.0f.\n",i+1,initial_points[i][0],initial_points[i][1]);
        std::cout << "alpha: " << alpha*180/M_PI << "; f(alpha): " << falpha << std::endl;
    }
    
    ///--------------get all solutions on [0,360)---------------
    std::cout << std::endl << "all solutions on [0,360):" << std::endl;
    for(double i = 0 ; i <360 ; i += 0.0001)
    {
        double t = i * M_PI /180;
        if(std::abs(f(t)) < 1e-5)
        {
            std::cout << i << " " << f(t) << std::endl;
            i = i + 1;
        }
    }

    return 0;
}