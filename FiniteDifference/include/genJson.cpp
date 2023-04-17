/**
 * @file genJson.cpp
 * @author MaCheng (Ch.Ma01@outlook.com)
 * @brief implement of genJson.h
 * @version 0.1
 * @date 2022-03-25
 * 
 * @copyright Copyright (c) 2022
 * 
 */


#include "genJson.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <queue>

genJson::genJson(double _h):h(_h), center(Vec(-1,-1)), radius(0)
{   
    if(!isInteger(1/h))
        std::cerr << "Invalid h." << std::endl;
    const int N = round(1/h);
    f_val.resize(N+1);
    for(int i = 0; i <= N ; i ++)
        f_val[i].resize(N+1);
};

bool genJson::isInDomain(Vec _p)
{
    if(dist(_p, center) <= radius)
        return false;
    if(_p[0] > 0 && _p[0] < 1 && _p[1] > 0 && _p[1] < 1)
        return true;
    return false;
}

bool genJson::isInSquare(Vec _p)
{
    if(_p[0]>=0 && _p[0] <= 1 && _p[1] >= 0 && _p[1] <= 1)
        return true;
    else
        return false;
}


int genJson::set_circ(Vec _center, double _radius)
{
    const int N = round(1/h);

    int visited[(N+1)*(N+1)] = {0};
    std::queue<int> BFSq;
    Vec a;
    for(int i = 0 ; i <= N ; i ++)
        for(int j = 0 ; j <= N ; j ++)
            if(dist(Vec(j*h, i*h), center) > _radius +  EPS)
            {
                BFSq.push(to1d(i,j));
                visited[to1d(i,j)] = 1;
                break;
            }
    if(BFSq.empty())
        return -1;
    while(!BFSq.empty())
    {
        auto t = BFSq.front();
        int i = to2d1(t);
        int j = to2d2(t);
        Vec p(j*h, i*h);
        if(i+1 <= N && !visited[to1d(i+1,j)])
        {
            visited[to1d(i+1,j)] = 1;
            if(dist(UP(p), _center) > _radius + EPS)
                BFSq.push(to1d(i+1, j));
        }
        if(j+1 <= N && !visited[to1d(i,j+1)])
        {
            visited[to1d(i,j+1)] = 1;
            if(dist(RIGHT(p), _center) > _radius + EPS)
                BFSq.push(to1d(i, j+1));
        }
        if(i-1 >= 0 && !visited[to1d(i-1,j)])
        {
            visited[to1d(i-1,j)] = 1;
            if(dist(DOWN(p), _center) > _radius + EPS)
                BFSq.push(to1d(i-1, j));
        }
        if(j-1 >= 0 && !visited[to1d(i,j-1)])
        {
            visited[to1d(i,j-1)] = 1;
            if(dist(LEFT(p), _center) > _radius + EPS)
                BFSq.push(to1d(i, j-1));
        }
        BFSq.pop();
    }

    int sum = 0;
    for(int i = 0 ; i <= N ; i ++)
        for(int j = 0 ; j <= N ; j ++)
            if(dist(Vec(j*h, i*h), _center) > _radius)
            {
                Vec p(i*h, j*h);
                if(dist(p, _center) > _radius + EPS && !visited[to1d(i,j)])
                    return -1;
                sum ++;
            }
    if( sum < 4)
        return -1;

    radius = _radius;
    center = _center;
    return 0;
}

Vec genJson::toGrid(Vec _p)
{
    Vec res;
    res = _p;
    if(_p[0] < 0 || _p[0] > 1 || _p[1] < 0 || _p[1] > 1)
        return _p;
    if(isInteger(_p[0]/h))
        _p[0] = round(_p[0]/h)*h;
    if(isInteger(_p[1]/h))
        _p[1] = round(_p[1]/h)*h;
    return _p;
}


std::vector<Vec> genJson::discrete_circ()
{
    std::vector<Vec> res;
    const int N = round(1/h);
    if(radius == 0)
        return res;
    for(int n = 0 ; n <= N; n += 1) 
    {
        double i = n * h;
        double t = std::abs(i - center[0]);
        if(t > radius + EPS) // seperate 
            continue;
        if(std::abs(t - radius) < EPS) // tangency 
            if(isInSquare(Vec(i, center[1])))
                res.push_back(toGrid(Vec(i, center[1])));
        t = std::sqrt(radius*radius - t*t); 
        if(isInSquare(Vec(i, center[1]+t)))
            res.push_back(toGrid(Vec(i, center[1]+t)));
        if(isInSquare(Vec(i, center[1]-t)))
            res.push_back(toGrid(Vec(i, center[1]-t)));
    }

    for(int n = 0 ; n <= N; n += 1) 
    {
        double i = n * h;
        double t = std::abs(i - center[1]);
        if(t > radius + EPS) // seperate 
            continue;
        if(std::abs(t - radius) < EPS) // tangency 
        {
            if(isInteger(center[0]/h))
                continue;
            if(isInSquare(Vec(center[0], i)))
                res.push_back(toGrid(Vec(center[0], i)));
            continue;
        }
        // intersect
        t = std::sqrt(radius*radius - t*t); 
        if(!isInteger((center[0]+t)/h) && isInSquare(Vec(center[0]+t, i)))
            res.push_back(toGrid(Vec(center[0]+t, i)));
        if(!isInteger((center[0]-t)/h) && isInSquare(Vec(center[0]-t, i)))
            res.push_back(toGrid(Vec(center[0]-t, i)));
    }
    // std::cout << res.size() << std::endl;
    return res;
}

void genJson::set_BV(Vec _p, Vec _n, double _val)
{
    auto it = BV.begin();
    for(; it != BV.end() ; it ++)
    {
        Vec p = it->p;
        if(std::abs(p[0]-_p[0]) < EPS && std::abs(p[1]-_p[1]) < EPS)    
            break;
    }
    if(it == BV.end())
        BV.push_back(BoundaryValues(_p, _n, _val));
    else
    {
        it -> n_vec = _n;
        it -> val = _val;
    }
    return;
}

void genJson::set_BV0(double (*u)(double, double))
{
    std::vector<Vec> set = discrete_circ();
    const int N = round(1/h);
    for(auto it = set.begin(); it != set.end(); it ++)
    {
        Vec p = *it;
        BV.push_back(BoundaryValues(p, Vec(0,0), u(p[0], p[1])));
    }
    for(int i = 0 ; i <= N ; i ++)
    {
        Vec p1(0, i*h), p2(1, i*h);
        if(dist(p1, center) > radius+EPS)
            BV.push_back(BoundaryValues(p1, Vec(0,0), u(p1[0], p1[1])));
        if(dist(p2, center) > radius+EPS)
            BV.push_back(BoundaryValues(p2, Vec(0,0), u(p2[0], p2[1])));
    }
    for(int j = 1 ; j <= N-1 ; j ++)
    {
        Vec p1(j*h, 0), p2(j*h, 1);
        if(dist(p1, center) > radius+EPS)
            BV.push_back(BoundaryValues(p1, Vec(0,0), u(p1[0], p1[1])));
        if(dist(p2, center) > radius+EPS)
            BV.push_back(BoundaryValues(p2, Vec(0,0), u(p2[0], p2[1])));
    }
    return;
}

void genJson::set_square_BV0(double (*u)(double, double))
{
    const int N = round(1/h);
    for(int i = 0 ; i <= N ; i ++)
    {
        Vec p1(0, i*h), p2(1, i*h);
        if(dist(p1, center) > radius+EPS)
            BV.push_back(BoundaryValues(p1, Vec(0,0), u(p1[0], p1[1])));
        if(dist(p2, center) > radius+EPS)
            BV.push_back(BoundaryValues(p2, Vec(0,0), u(p2[0], p2[1])));
    }
    for(int j = 1 ; j <= N-1 ; j ++)
    {
        Vec p1(j*h, 0), p2(j*h, 1);
        if(dist(p1, center) > radius+EPS)
            BV.push_back(BoundaryValues(p1, Vec(0,0), u(p1[0], p1[1])));
        if(dist(p2, center) > radius+EPS)
            BV.push_back(BoundaryValues(p2, Vec(0,0), u(p2[0], p2[1])));
    }
    return;
}

void genJson::set_circle_BV0(double (*u)(double, double))
{
    std::vector<Vec> set = discrete_circ();
    const int N = round(1/h);
    for(auto it = set.begin(); it != set.end(); it ++)
    {
        Vec p = *it;
        BV.push_back(BoundaryValues(p, Vec(0,0), u(p[0], p[1])));
    }
    return;
}

void genJson::set_BV1(double (*ux)(double, double), double (*uy)(double, double), Vec _n)
{
    std::vector<Vec> set = discrete_circ();
    const int N = round(1/h);
    if(_n == Vec(0,0))  
        std::cerr << "Invalid n_vec" << std::endl;
    _n = _n/_n.norm();
    for(auto it = set.begin(); it != set.end(); it ++)
    {
        Vec p = *it;
        double val = ux(p[0],p[1])*_n[0] + uy(p[0],p[1])*_n[1];
        BV.push_back(BoundaryValues(p ,_n , val));
    }
    for(int i = 0 ; i <= N ; i ++)
    {
        Vec p1(0, i*h), p2(1, i*h);
        if(dist(p1, center) > radius+EPS)
            BV.push_back(BoundaryValues(p1, _n, ux(p1[0],p1[1])*_n[0]+uy(p1[0],p1[1])*_n[1]));
        if(dist(p2, center) > radius+EPS)
            BV.push_back(BoundaryValues(p2, _n, ux(p2[0],p2[1])*_n[0]+uy(p2[0],p2[1])*_n[1]));
    }
    for(int j = 1 ; j <= N-1 ; j ++)
    {
        Vec p1(j*h, 0), p2(j*h, 1);
        if(dist(p1, center) > radius+EPS)
            BV.push_back(BoundaryValues(p1, _n, ux(p1[0],p1[1])*_n[0]+uy(p1[0],p1[1])*_n[1]));
        if(dist(p2, center) > radius+EPS)
            BV.push_back(BoundaryValues(p2, _n, ux(p2[0],p2[1])*_n[0]+uy(p2[0],p2[1])*_n[1]));
    }
    return;
}

void genJson::set_circle_BV1(double (*ux)(double, double), double (*uy)(double, double), Vec _n)
{
    std::vector<Vec> set = discrete_circ();
    const int N = round(1/h);
    if(_n == Vec(0,0))  
        std::cerr << "Invalid n_vec" << std::endl;
    _n = _n/_n.norm();
    for(auto it = set.begin(); it != set.end(); it ++)
    {
        Vec p = *it;
        double val = ux(p[0],p[1])*_n[0] + uy(p[0],p[1])*_n[1];
        BV.push_back(BoundaryValues(p ,_n , val));
    }
    return;
}

void genJson::set_square_BV1(double (*ux)(double, double), double (*uy)(double, double), Vec _n)
{
    std::vector<Vec> set = discrete_circ();
    const int N = round(1/h);
    if(_n == Vec(0,0))  
        std::cerr << "Invalid n_vec" << std::endl;
    _n = _n/_n.norm();
    for(int i = 0 ; i <= N ; i ++)
    {
        Vec p1(0, i*h), p2(1, i*h);
        if(dist(p1, center) > radius+EPS)
            BV.push_back(BoundaryValues(p1, _n, ux(p1[0],p1[1])*_n[0]+uy(p1[0],p1[1])*_n[1]));
        if(dist(p2, center) > radius+EPS)
            BV.push_back(BoundaryValues(p2, _n, ux(p2[0],p2[1])*_n[0]+uy(p2[0],p2[1])*_n[1]));
    }
    for(int j = 1 ; j <= N-1 ; j ++)
    {
        Vec p1(j*h, 0), p2(j*h, 1);
        if(dist(p1, center) > radius+EPS)
            BV.push_back(BoundaryValues(p1, _n, ux(p1[0],p1[1])*_n[0]+uy(p1[0],p1[1])*_n[1]));
        if(dist(p2, center) > radius+EPS)
            BV.push_back(BoundaryValues(p2, _n, ux(p2[0],p2[1])*_n[0]+uy(p2[0],p2[1])*_n[1]));
    }
    return;
}


int genJson::set_fval(double (*f)(double, double))
{
    const int N = round(1/h);
    for(int i = 0 ; i <= N ; i ++)
    {
        for(int j = 0 ; j <= N ; j ++)
        {
            if(!isInDomain(Vec(j*h, i*h)))
            {
                f_val[i][j] = D_MAX;
                continue;
            }
            f_val[i][j] = f(j*h, i*h);
        }
    }
    return 0;
}

int genJson::write_json(std::string _path)
{
    Json::Value root;
    Json::StyledWriter writer;
    std::ofstream outjson;

    const int N = round(1/h);

    root["h"] = h; // grid size

    // circle
    Json::Value json_circ;
    json_circ["center"].append(center[0]);
    json_circ["center"].append(center[1]);
    json_circ["radius"] = radius;
    root["circle"] = json_circ;

    // value of f
    Json::Value json_f_val;
    for(int i = 0; i <= N ; i ++ )
        for(int j = 0 ; j <= N ; j ++)
        {
            Json::Value json_tmp;
            json_tmp.append(i);
            json_tmp.append(j);
            json_tmp.append(f_val[i][j]);
            json_f_val.append(json_tmp);
        }
    root["f_val"] = json_f_val;

    // BV
    Json::Value json_BV;
    for(auto it = BV.begin(); it != BV.end() ; it ++)
    {
        Json::Value json_tmp;
        auto t = *it;
        json_tmp[0].append(t.p[0]);
        json_tmp[0].append(t.p[1]);
        json_tmp[1].append(t.n_vec[0]);
        json_tmp[1].append(t.n_vec[1]);
        json_tmp[2] = t.val;
        json_BV.append(json_tmp);
    }
    root["BV"] = json_BV;

    outjson.open(_path);
    outjson << writer.write(root);
    outjson.close();
    return 0;
}

///////////////////test
/*
double f(double x, double y)
{
    return std::exp(y + std::sin(x))*std::sin(x) - std::exp(y + std::sin(x)) - std::exp(y + std::sin(x))*std::cos(x)*std::cos(x);
    // return 0;
}

double u(double x, double y)
{
    return std::exp(y + std::sin(x));
    // return 1;
}

int main()
{
    double const h = 0.1;
    int const N = 1/h;
    std::string path = "./first.json";
    genJson mycond(h);

    mycond.set_circ(Vec(0.5,0.5), 0.15);
    mycond.set_fval(f);
    mycond.set_BV0(u);
    mycond.write_json(path);
    // mycond.set_BV_square(0, u);
    // mycond.set_BV_circ(u);
    return 0;
}*/
