/**
 * @file global.cpp
 * @author MaCheng (Ch.Ma01@outlook.com)
 * @brief implement of global.h
 * @version 0.1
 * @date 2022-03-25
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include "global.h"

bool isInteger(double t)
{
    if(std::abs(int(t + 0.5) - t) < EPS)
        return true;
    return false;
}

bool isInSquare(Vec _p)
{
    if(_p[0] >= 0-EPS && _p[0] <= 1+EPS && _p[1] >= 0-EPS && _p[1] <= 1+EPS)
        return true;
    return false;
}

double round(double t){ return int(t+0.5); };

double max(double x, double y)
{
    if(x > y) 
        return x;
    return y;
};