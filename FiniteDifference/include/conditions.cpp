/**
 * @file conditions.cpp
 * @author MaCheng (Ch.Ma01@outlook.com)
 * @brief implement of conditions.h
 * @version 0.1
 * @date 2022-03-03
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <assert.h>
#include "conditions.h" 
#include "Vec2d.h"

/**
 * @brief Determine whether a float number is an integer.
 * 
 * @param t given float number.
 * @return true: t is an integer.
 */

Conditions::Conditions(std::string _path)
{
    std::ifstream infile(_path);
    assert(infile.good());
    Json::Reader reader;
    Json::Value root;
    if(!reader.parse(infile, root))
    {
        std::cerr << "read error!" << std::endl;
        return;
    }
    
    h = root["h"].asDouble(); // h

    // setting the size of f_val
    const int N = round(1/h);
    f_val.resize(N+1);
    for(auto &t: f_val)
        t.resize(N+1);

    // circle
    Json::Value circle = root["circle"];
    center[0] = circle["center"][0].asDouble();
    center[1] = circle["center"][1].asDouble();
    radius = circle["radius"].asDouble();

    // f_val
    Json::Value json_fval = root["f_val"];
    for(auto t: json_fval)
        f_val[t[0].asInt()][t[1].asInt()] = t[2].asDouble();

    // BV
    Json::Value json_BV = root["BV"];
    for(auto t: json_BV)   
    {
        Vec p(t[0][0].asDouble(), t[0][1].asDouble());
        Vec n(t[1][0].asDouble(), t[1][1].asDouble());
        BV.push_back(BoundaryValues(p, n, t[2].asDouble()));
    }
}

struct Conditions::BoundaryValues &Conditions::find_BV(Vec _p)
{
    for(auto it = BV.begin(); it != BV.end() ; it ++)
    {
        Vec p = it->p;
        if(std::abs(p[0]-_p[0]) < EPS && std::abs(p[1]-_p[1]) < EPS)    
            return *it;   
    }
}

int Conditions::find_i(int i, double x)
{
    int j = int(x/h);
    const int N = round(1/h);
    for(int k = 1 ; k <= N ; k ++)
    {
        Vec u1(j*h, (i+k)*h), u2((j+1)*h, (i+k)*h), u3(x, (i+k)*h);
        Vec d1(j*h, (i-k)*h), d2((j+1)*h, (i-k)*h), d3(x, (i-k)*h);
        if(i+k <= N && dist(u1, center) > radius-EPS && dist(u2, center) > radius-EPS && dist(u3, center) > radius-EPS)
            return i+k;
        if(i-k >= 0 && dist(d1, center) > radius-EPS && dist(d2, center) > radius-EPS && dist(d3, center) > radius-EPS)
            return i-k;
    }     
    return -1;  
}

int Conditions::find_j(int j, double y)
{
    int i = int(y/h);
    const int N = round(1/h);
    for(int k = 1 ; k <= N ; k ++)
    {
        Vec l1((j-k)*h, i*h), l2((j-k)*h, (i+1)*h), l3((j-k)*h, y);
        Vec r1((j+k)*h, i*h), r2((j+k)*h, (i+1)*h), r3((j+k)*h, y);
        if(j+k <= N && dist(r1, center) > radius-EPS && dist(r2, center) > radius-EPS && dist(r3, center) > radius-EPS)
            return j+k;
        if(j-k >= 0 && dist(l1, center) > radius-EPS && dist(l2, center) > radius-EPS && dist(l3, center) > radius-EPS)
            return j-k;
    }     
    return -1;  
}