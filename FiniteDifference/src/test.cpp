/**
 * @file test.cpp
 * @author OuyangShangke
 * @brief to execute all the tests
 * @version 0.1
 * @date 2022-03-25
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include "FD.h"
#include "conditions.h"
#include "global.h"
#include "genJson.h"
#include <iostream>
#include <math.h>

typedef double (*FUNCTION) (double x, double y);

double u1(double x, double y)
{
    return std::exp(y + std::sin(x));
}

int start_test(int index, FUNCTION u)
{
    double n_set[4] = {8, 16, 32, 64};
    double1d n_list(n_set, n_set+4);
    std::cout << "test index: " << index << std::endl;
    double1d inf,one,two;
    for(auto n: n_list)
    {
        std::cout << "N = " << int(n) << std::endl;
        cstring path = cstring("../json/f")+std::to_string(index)+cstring("_N")+std::to_string(int(n))+cstring(".json");
        cstring outpath = cstring("../plot/f")+std::to_string(index)+cstring("_N")+std::to_string(int(n))+cstring(".m");
        Conditions mycond(path);
        FDsolver mysolver(mycond);
        mysolver.general_solver();
        inf.push_back(mysolver.get_err(u, 0));
        one.push_back(mysolver.get_err(u, 1));
        two.push_back(mysolver.get_err(u, 2));
        double r_inf = 0, r1 = 0, r2 = 0;
        if(n != 8)
        {
            r_inf = -(std::log(*(inf.end()-1)) - std::log(*(inf.end()-2)))/std::log(2);
            r1 = -(std::log(*(one.end()-1)) - std::log(*(one.end()-2)))/std::log(2);
            r2 = -(std::log(*(two.end()-1)) - std::log(*(two.end()-2)))/std::log(2);
        }
        std::cout <<"inf-norm err: " << *(inf.end()-1) << "; convergence rate: " << r_inf << std::endl;
        std::cout <<"1-norm err: " << *(one.end()-1) << "; convergence rate: " << r1 << std::endl;
        std::cout <<"2-norm err: " << *(two.end()-1) << "; convergence rate: " << r2 << std::endl;
        mysolver.get_err(u, outpath, index);
    }
}

int main()
{
    FUNCTION pFuncList[] = {&u1};
    int test_num = 9;
    for(int i = 1 ; i <= test_num ; i ++)
        start_test(i, u1);
        // start_test(i, pFuncList[i-1]);
    // Conditions mycond(Dirichlet, regular, h, f);
    return 0;
}