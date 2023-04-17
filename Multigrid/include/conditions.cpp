/**
 * @file conditions.cpp
 * @author OuyangShangke
 * @brief implement of conditions.h
 * @version 0.1
 * @date 2022-03-03
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include "conditions.h" 

bool isPower2(int _t)
{
    if(_t <= 0)
        return false;
    while(_t % 2 == 0)
        _t /= 2;
    if(_t == 1)
        return true;
    return false;
}

Conditions<1>::Conditions(std::string _path)
{
    std::ifstream infile(_path);
    assert(infile.good());
    Json::Value root;

    Json::CharReaderBuilder builder;
	builder["collectComments"] = true;
	Json::String errs;
	assert(parseFromStream(builder, infile, &root, &errs) && "error when parse json");

    // if(!reader.parse(infile, root))
    // {
    //     std::cerr << "read error!" << std::endl;
    //     return;
    // }

    N = root["N"].asInt(); // h
    assert(N >= 8 && isPower2(N));

    // setting the size of f_val
    f_val.resize(N+1);

    // f_val
    Json::Value json_fval = root["f_val"];
    for(auto t: json_fval)
        f_val[t[0].asInt()] = t[1].asDouble();

    // BV
    Json::Value json_BV = root["BV"];
    for(auto t: json_BV)   
    {
        BV[t[0].asInt()][0] = t[1].asDouble(); // value
        BV[t[0].asInt()][1] = t[2].asInt(); // type of BV
    }
}

Conditions<2>::Conditions(std::string _path)
{
    std::ifstream infile(_path);
    assert(infile.good());
    Json::Value root;
    
    Json::CharReaderBuilder builder;
	builder["collectComments"] = true;
	Json::String errs;
	assert(parseFromStream(builder, infile, &root, &errs) && "error when parse json");
    
    N = root["N"].asDouble(); // h
    assert(N >= 8 && isPower2(N) && "invalid spectrum.");

    // setting the size of f_val
    f_val.resize(N+1);
    for(auto &t: f_val)
        t.resize(N+1);

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

    // if it's regular domain
    isRegular = root["isRegular"].asInt();
    if(!isRegular) // irregular
    {
        loc.resize(N+1);
        for(auto &t: loc)
            t.resize(N+1);
        Json::Value json_loc = root["loc"];
        for(auto t: json_loc)
            loc[t[0].asInt()][t[1].asInt()] = t[2].asInt();
    }
    
}

const Conditions<2>::BoundaryValues& Conditions<2>::find_BV(Vec _p)
{
    for(auto it = BV.begin(); it != BV.end() ; it ++)
    {
        Vec p = it->p;
        if(std::abs(p[0]-_p[0]) < EPS && std::abs(p[1]-_p[1]) < EPS)    
            return *it;   
    }
    assert(0 && "bug");
}

/**
 * @brief Whether point p3 is between p1 and p2
 * 
 * @param _p1 
 * @param _p2 
 * @param _p3 
 * @return true between
 * @return false
 */
bool isBetween(Vec _p1, Vec _p2, Vec _p3)
{
    if(isEqual(_p1[0], _p2[0])) // same x
    {
        if(!isEqual(_p1[0], _p3[0]))
            return false;
        double a = _p1[1], b = _p2[1];
        if(a > b)
            std::swap(a,b);
        if(a < _p3[1] && _p3[1] < b) 
            return true;
        return false;
    }
    if(isEqual(_p1[1], _p2[1])) // same y
    {
        if(!isEqual(_p1[1], _p3[1]))
            return false;
        double a = _p1[0], b = _p2[0];
        if(a > b)
            std::swap(a,b);
        if(a < _p3[0] && _p3[0] < b) 
            return true;
        return false;
    }
    return false;
}

const Conditions<2>::BoundaryValues& Conditions<2>::find_BV(Vec _p1, Vec _p2)
{
    int res = -1;
    double min = 100;
    for(int i = 0 ; i < BV.size() ; i ++)
    {
        Vec p = BV[i].p;
        if(isBetween(_p1, _p2, p))    
        {
            double t = dist(p, _p1);
            if(t < min)
            {
                min = t;
                res = i;
            }
        }
    }
    // std::cout << _p1 << " " << _p2 << std::endl;
    assert(res != -1);
    return BV[res];
}
