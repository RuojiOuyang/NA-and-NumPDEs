#include "EquationSolver.cpp"

double f1(double x)
{
    return sin(x / 2) - 1;
}

double f2(double x)
{
    return exp(x) - tan(x);
}

double f3(double x)
{
    return pow(x , 3) - 12 * pow(x , 2) + 3 * x + 1;
}

int main(int argc , char* argv[])
{
    double x , fx;
    ///(1)
    double initial_values1[2][2] = {{0 , M_PI / 2},{4*M_PI, 5*M_PI}};
    class secant solver1(f1 , 0 , 0);
    std::cout << "(1)" << std::endl;
    for(int i = 0 ; i < 2 ;  i ++)
    {   
        double x0 = initial_values1[i][0];
        double x1 = initial_values1[i][1];
        solver1.set_initial(x0, x1);
        printf("%d. initial points: %.5f, %.5f.\n",i+1,x0,x1);
        std::cout << "solution: " << solver1.solve() << std::endl;
    }

    ///(2)
    double initial_values2[2][2] = {{1 , 1.4},{-2 , -3}};
    class secant solver2(f2 , 0 , 0);
    std::cout << "(2)" << std::endl;
    for(int i = 0 ; i < 2 ;  i ++)
    {   
        double x0 = initial_values2[i][0];
        double x1 = initial_values2[i][1];
        solver2.set_initial(x0, x1);
        printf("%d. initial points: %.5f, %.5f.\n",i+1,x0,x1);
        std::cout << "solution: " << solver2.solve() << std::endl;
    }

    ///(3)
    double initial_values3[2][2] = {{ 0 , -0.5},{0 , 0.5}};
    class secant solver3(f3 , 0 , 0);
    std::cout << "(3)" << std::endl;
    for(int i = 0 ; i < 2 ;  i ++)
    {   
        double x0 = initial_values3[i][0];
        double x1 = initial_values3[i][1];
        solver3.set_initial(x0, x1);
        printf("%d. initial points: %.5f, %.5f.\n",i+1,x0,x1);
        std::cout << "solution: " << solver3.solve() << std::endl;
    }
    return 0;
}
