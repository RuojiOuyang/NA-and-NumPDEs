/**
 * @file global.cpp
 * @author OuyangShangke
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

bool isEqual(double a, double b)
{
    return (std::abs(a-b)<EPS);
}

int to1d(int i , int  j , int N)
{
    return i*(N+1)+j;
}

double round(double t){ return int(t+0.5); };

double max(double x, double y)
{
    if(x > y) 
        return x;
    return y;
};

