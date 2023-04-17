/**
 * @file genJson.cpp
 * @author OuyangShangke
 * @brief implement of genJson.h
 * @version 0.1
 * @date 2022-03-25
 * 
 * @copyright Copyright (c) 2022
 * 
 */


#include "genJson.h"

/**
 * 1-dimension
 */
int genJson<1>::set_fval(FUNC1d _f)
{
    double h = 1/double(N);
    f_val.resize(N+1);
    for(int i = 1 ; i <= N-1 ; i ++)
        f_val[i] = _f(i*h);
    return 0;
}

int genJson<1>::set_BV(int _p, int _type, double _val)
{
    BV[_p][0] = _val;
    BV[_p][1] = _type;
    return 0;
}

int genJson<1>::write_json(std::string _path)
{
    Json::Value root;
    Json::StreamWriterBuilder writerBuilder;
    std::ofstream outjson;

    root["dimension"] = 1;

    root["N"] = N; // grid size

    // value of f
    Json::Value json_f_val;
    for(int i = 0; i <= N ; i ++ )
    {
        Json::Value json_tmp;
        json_tmp.append(i);
        json_tmp.append(f_val[i]);
        json_f_val.append(json_tmp);
    }
    root["f_val"] = json_f_val;

    // BV
    Json::Value json_BV;
    for(int i = 0 ; i <= 1 ; i ++)
    {
        Json::Value json_tmp;
        json_tmp.append(i);
        json_tmp.append(BV[i][0]);
        json_tmp.append(BV[i][1]);
        json_BV.append(json_tmp);
    }
    root["BV"] = json_BV;

    outjson.open(_path);
    std::unique_ptr<Json::StreamWriter> jsonWriter(writerBuilder.newStreamWriter());
    jsonWriter->write(root, &outjson);
    outjson.close();
    return 0;
}

/**
 * 2-dimension
 */

// genJson<2>::genJson(int _N, bool _isReg):N(_N), isRegular(_isReg)
// {
//     f_val.resize(N+1);
//     for(auto &t: f_val)
//         t.resize(N+1);
//     if(!isRegular)
//     {
//         loc.resize(N+1);
//         for(auto &t: loc)
//             t.resize(N+1);
//     }
// }

int genJson<2>::set_sine_domain(double _spec, FUNC2d _f, FUNC2d _u)
{
    assert(0<_spec && _spec < 1 && "invalid spectrum");
    isRegular = false;
    f_val.resize(N+1);
    for(auto &t: f_val)
        t.resize(N+1);
    loc.resize(N+1);
    for(auto &t: loc)
        t.resize(N+1);
    
    double h = 1/double(N);
    // set f_val
    for(int i = 0 ; i <= N ; i ++)
        for(int j = 0 ; j <= N ; j ++)
        {
            double x = j*h, y = i*h;
            if(y < _spec*std::sin(M_PI*x) - EPS) // out domain
            {
                loc[i][j] = 0;
            }
            else if(i == 0 || j == 0 || i == N || j == N || isEqual(_spec*std::sin(M_PI*x), y)) // on boundary
            {
                loc[i][j] = 1;
                BV.push_back(BoundaryValues(Vec(x,y), Vec(0,0), _u(x,y)));
            }
            else // in domain
            {
                loc[i][j] = 2;
                f_val[i][j] = _f(x,y);
            }
        }
    // set boundary on sine curve
    for(int j = 1 ; j < N ; j ++)
    {
        double x = j*h, y = _spec*std::sin(M_PI*x);
        if(isInteger(y/h))
            continue;
        BV.push_back(BoundaryValues(Vec(x,y), Vec(0,0), _u(x,y)));
    }    
    for(int i = 1 ; i <= N ; i ++)
    {
        double y = i*h;
        if(y > _spec - EPS) // y >= _spec
            break;
        else // y < _spec
        {
            double t = y/_spec;
            t = std::asin(t);
            double x1 = t/M_PI, x2 = 1 - x1;
            if(!isInteger(x1/h))
            {
                BV.push_back(BoundaryValues(Vec(x1,y), Vec(0,0), _u(x1,y)));
                BV.push_back(BoundaryValues(Vec(x2,y), Vec(0,0), _u(x2,y)));
            }
        }
    }   
    return 0;
}

int genJson<2>::set_fval(FUNC2d _f)
{
    double h = 1/double(N);
    f_val.resize(N+1);
    for(auto &t: f_val)
        t.resize(N+1);
    for(int i = 0 ; i <= N ; i ++)
        for(int j = 0 ; j <= N ; j ++)
            f_val[i][j] = _f(j*h, i*h);
    return 0;
}

int genJson<2>::set_BV0(FUNC2d _u)
{
    double h = 1/double(N);
    for(int i = 0 ; i <= N ; i ++)
    {
        Vec p1(0, i*h), p2(1, i*h);
        BV.push_back(BoundaryValues(p1, Vec(0,0), _u(p1[0], p1[1])));
        BV.push_back(BoundaryValues(p2, Vec(0,0), _u(p2[0], p2[1])));
    }
    for(int j = 1 ; j <= N-1 ; j ++)
    {
        Vec p1(j*h, 0), p2(j*h, 1);
        BV.push_back(BoundaryValues(p1, Vec(0,0), _u(p1[0], p1[1])));
        BV.push_back(BoundaryValues(p2, Vec(0,0), _u(p2[0], p2[1])));
    }
    return 0;
}

int genJson<2>::set_BV1(FUNC2d _ux, FUNC2d _uy, Vec _n)
{
    double h = 1/double(N);
    assert(!(_n == Vec(0,0)) && "invalid n_vec");  
    _n = _n/_n.norm();
    for(int i = 0 ; i <= N ; i ++)
    {
        Vec p1(0, i*h), p2(1, i*h);
        BV.push_back(BoundaryValues(p1, _n, _ux(p1[0],p1[1])*_n[0]+_uy(p1[0],p1[1])*_n[1]));
        BV.push_back(BoundaryValues(p2, _n, _ux(p2[0],p2[1])*_n[0]+_uy(p2[0],p2[1])*_n[1]));
    }
    for(int j = 1 ; j <= N-1 ; j ++)
    {
        Vec p1(j*h, 0), p2(j*h, 1);
        BV.push_back(BoundaryValues(p1, _n, _ux(p1[0],p1[1])*_n[0]+_uy(p1[0],p1[1])*_n[1]));
        BV.push_back(BoundaryValues(p2, _n, _ux(p2[0],p2[1])*_n[0]+_uy(p2[0],p2[1])*_n[1]));
    }
    return 0;
}


int genJson<2>::set_BV(Vec _p, Vec _n, double _val)
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
    return 0;
}

int genJson<2>::write_json(std::string _path)
{
    Json::Value root;
    Json::StreamWriterBuilder writerBuilder;
    std::ofstream outjson;

    root["dimension"] = 2;

    root["isRegular"] = int(isRegular);

    root["N"] = N; // the number of cells

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
    //loc
    if(!isRegular)
    {
        Json::Value json_loc;
        for(int i = 0; i <= N ; i ++ )
            for(int j = 0 ; j <= N ; j ++)
            {
                Json::Value json_tmp;
                json_tmp.append(i);
                json_tmp.append(j);
                json_tmp.append(loc[i][j]);
                json_loc.append(json_tmp);
            }
        root["loc"] = json_loc;
    }

    outjson.open(_path);
    std::unique_ptr<Json::StreamWriter> jsonWriter(writerBuilder.newStreamWriter());
    jsonWriter->write(root, &outjson);
    outjson.close();
    return 0;
}

// template class genJson<1>;
// template class genJson<2>;

////////test
// double f(double x, double y)
// {return x+y;}
// double u(double x , double y )
// {return x - y;}

// int main()
// {
//     genJson<2> my(8, false);
//     my.set_sine_domain(1/16.0, f, u);
//     my.write_json(string("./my.json"));
//     return 0;
// }

