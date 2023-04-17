/**
 * @file myJson6.cpp
 * @author OuyangShangke
 * @brief generate json for: 
 * u(x,y) = exp(y + sin(x)) \n
 * domain: square minus a circle whoes center is at the center of square \n 
 * boundary conditions: Mixed. \n
 * - boundary of circle and corners of square give Dirichlet conditions
 * @version 0.1
 * @date 2022-03-23
 * 
 * @copyright Copyright (c) 2022
 * 
 */


#include "FD.h"
#include "conditions.h"
#include "global.h"
#include "genJson.h"
#include <iostream>
#include <jsoncpp/json/json.h>
#include <math.h>


double f(double x, double y)
{
    return std::exp(y + std::sin(x))*std::sin(x) - std::exp(y + std::sin(x)) - std::exp(y + std::sin(x))*std::cos(x)*std::cos(x);
}

double u(double x, double y)
{
    return std::exp(y + std::sin(x));
}

double ux(double x, double y)
{
    return std::exp(y + std::sin(x))*std::cos(x);
}

double uy(double x, double y)
{
    return std::exp(y + std::sin(x));
}

int main()
{
    double n_set[4] = {8, 16, 32, 64};
    double1d n_list(n_set, n_set+4);
    
    for(auto n: n_list)
    {
        genJson myjson(1/n);
        cstring path = cstring("../json/f6_N")+std::to_string(int(n))+cstring(".json");
        myjson.set_circ(Vec(0.5,0.5), 0.15);
        myjson.set_fval(f);
        myjson.set_circle_BV0(u);
        myjson.set_square_BV1(ux, uy, Vec(1,2));
        myjson.set_BV(Vec(0, 0), Vec(0,0), u(0, 0));
        myjson.set_BV(Vec(0, 1), Vec(0,0), u(0, 1));
        myjson.set_BV(Vec(1, 0), Vec(0,0), u(1, 0));
        myjson.set_BV(Vec(1, 1), Vec(0,0), u(1, 1));
        myjson.write_json(path);
    }
    return 0;
}