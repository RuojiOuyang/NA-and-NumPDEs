/**
 * @file operator.cpp
 * @author OuyangShangke
 * @brief implement of operator.h
 * @version 0.1
 * @date 2022-04-14
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include "operator.h"

// OP<1>
VectorXd OP<1>::restriction(restriction_type _type, const VectorXd& _vec)
{
    if(_type == full_weighting)
    {
        const int N = _vec.size() - 1;
        assert(N % 2 == 0);
        VectorXd res(N/2 + 1);
        res(0) = _vec(0);
        res(N/2) = _vec(N);
        for(int j = 1 ; j < N/2 ; j ++)
            res[j] = (_vec[2*j-1]+2.0*_vec[2*j]+_vec[2*j+1])/4.0;
        return res;
    }
    else if(_type == injection)
    {
        const int N = _vec.size() - 1;
        // std::cout << std::endl << N << std::endl;
        assert(N % 2 == 0);
        VectorXd res(N/2 + 1);
        for(int j = 0 ; j <= N/2 ; j ++)
            res[j] = _vec[2*j];
        return res;
    }   
    assert(0);
}

VectorXd OP<1>::interpolation(interpolation_type _type, const VectorXd& _vec)
{
    const int N = (_vec.size()-1)*2;
    VectorXd res(N+1);
    if(_type == linear)
    {
        for(int j = 0 ; j <= N ; j ++)
        {
            if(j%2 == 0)
                res[j] = _vec[j/2];
            else
                res[j] = (_vec[(j-1)/2]+_vec[(j+1)/2])/2.0;
        }
        return res;
    }
    else if(_type == quadratic)
    {
        for(int j = 0 ; j <= N ; j ++)
        {
            if(j%2 == 0)
                res[j] = _vec[j/2];
            else
            {
                double y1, y2, y3;
                if(j != N - 1)
                { 
                    int t = (j-1)/2;
                    y1 = _vec(t);
                    y2 = _vec(t+1);
                    y3 = _vec(t+2);
                }
                else
                {
                    int t = (j+1)/2;
                    y1 = _vec(t);
                    y2 = _vec(t-1);
                    y3 = _vec(t-2);
                }
                res[j] = 3/8.0*y1+3/4.0*y2-1/8.0*y3;
            }
        }
        return res;
    }
    assert(0);
}

// OP<2>

// OP<2>::OP(int _N):N_(_N)
// {
//     // uint1d tmp(N+1, 1);
//     // for(int i = 0 ; i <= N ; i ++)
//     //     loc.push_back(tmp);
//     ;
// }


VectorXd OP<2>::restriction(restriction_type _type, const VectorXd& _vec)
{
    if(_type == full_weighting)
    {
        int t = _vec.size();
        assert(isInteger(std::sqrt(t)));
        const int N = round(std::sqrt(t)) - 1;
        assert(N % 2 == 0);
        VectorXd res(SQUARE(N/2 + 1));
        

        for(int i = 0 ; i <= N/2 ; i ++)
            for(int j = 0 ; j <= N/2 ; j ++)
            {
                int s = 0;
                double sum = 0;
                if(i != 0)
                {
                    sum += _vec(to1d(i*2,j*2,N))+_vec(to1d(i*2-1,j*2,N));
                    s ++;
                }
                if(i != N/2)
                {
                    sum += _vec(to1d(i*2,j*2,N))+_vec(to1d(i*2+1,j*2,N));
                    s ++;
                }
                if(j != 0)
                {
                    sum += _vec(to1d(i*2,j*2,N))+_vec(to1d(i*2,j*2-1,N));
                    s ++;
                }
                if(j != N/2)
                {
                    sum += _vec(to1d(i*2,j*2,N))+_vec(to1d(i*2,j*2+1,N));
                    s ++;
                }
                res(to1d(i,j,N/2)) = sum/double(s*2);
            }
        return res;
    }
    else if(_type == injection)
    {
        int t = _vec.size();
        assert(isInteger(std::sqrt(t)));
        const int N = round(std::sqrt(t)) - 1;
        assert(N % 2 == 0);
        VectorXd res(SQUARE(N/2 + 1));
        for(int i = 0 ; i <= N/2 ; i ++)
            for(int j = 0 ; j <= N/2 ; j ++)
                res[to1d(i,j,N/2)] = _vec[to1d(i*2,j*2,N)];
        return res;
    }
    std::cout << _type << std::endl;
    assert(0);
}

VectorXd OP<2>::interpolation(interpolation_type _type, const VectorXd& _vec)
{
    int t = _vec.size();
    VectorXd vec = _vec;
    assert(isInteger(std::sqrt(t)));
    const int originN = cond.N;
    const int N = (round(std::sqrt(t)) - 1)*2;
    VectorXd res(SQUARE(N+1));
    /*
    if(!cond.isRegular) // irregular
    {
        // std::cout << "adjust" << std::endl;
        int n = N/2;
        double h = 1/double(n);
        int s = originN/n;
        for(int i = 0 ; i <= n ; i ++)
            for(int j = 0 ; j <= n ; j ++)
            {
                if(cond.loc[i*s][j*s] == 0) //outside
                {
                    // if(j+1 <= n && cond.loc[i*s][(j+1)*s] == 2)
                    // {
                    //     auto bv = cond.find_BV(Vec(j*h, i*h), Vec((j+1)*h, i*h));
                    //     double alpha = (bv.p[0] - j*h)/((j+1)*h - bv.p[0]);
                    //     vec[to1d(i,j,n)] = -alpha*vec[to1d(i,j+1,n)];
                    // }
                    // else if(j-1 >= 0 && cond.loc[i*s][(j-1)*s] == 2)
                    // {
                    //     auto bv = cond.find_BV(Vec(j*h, i*h), Vec((j-1)*h, i*h));
                    //     double alpha = (j*h - bv.p[0])/(bv.p[0] - (j-1)*h);
                    //     vec[to1d(i,j,n)] = -alpha*vec[to1d(i,j-1,n)];
                    // }
                    // else if(i+1 <= n && cond.loc[(i+1)*s][j*s] == 2)
                    // {
                    //     auto bv = cond.find_BV(Vec(j*h, i*h), Vec(j*h, (i+1)*h));
                    //     double alpha = (bv.p[1] - i*h)/((i+1)*h - bv.p[1]);
                    //     vec[to1d(i,j,n)] = -alpha*vec[to1d(i+1,j,n)];
                    // }
                    // else if(i-1 >= 0 && cond.loc[(i-1)*s][j*s] == 2)
                    // {
                    //     auto bv = cond.find_BV(Vec(j*h, i*h), Vec(j*h, (i-1)*h));
                    //     double alpha = (i*h - bv.p[1])/(bv.p[1] - (i-1)*h);
                    //     vec[to1d(i,j,n)] = -alpha*vec[to1d(i-1,j,n)];
                    // }
                    int tmp = 0;
                    double sum = 0;
                    if(j+1 <= n && cond.loc[i*s][(j+1)*s] == 2)
                    {
                        auto bv = cond.find_BV(Vec(j*h, i*h), Vec((j+1)*h, i*h));
                        double alpha = (bv.p[0] - j*h)/((j+1)*h - bv.p[0]);
                        sum += -alpha*vec[to1d(i,j+1,n)];
                        tmp++;
                    }
                    if(j-1 >= 0 && cond.loc[i*s][(j-1)*s] == 2)
                    {
                        auto bv = cond.find_BV(Vec(j*h, i*h), Vec((j-1)*h, i*h));
                        double alpha = (j*h - bv.p[0])/(bv.p[0] - (j-1)*h);
                        sum += -alpha*vec[to1d(i,j-1,n)];
                        tmp++;
                    }
                    if(i+1 <= n && cond.loc[(i+1)*s][j*s] == 2)
                    {
                        auto bv = cond.find_BV(Vec(j*h, i*h), Vec(j*h, (i+1)*h));
                        double alpha = (bv.p[1] - i*h)/((i+1)*h - bv.p[1]);
                        sum += -alpha*vec[to1d(i+1,j,n)];
                        tmp++;
                    }
                    if(i-1 >= 0 && cond.loc[(i-1)*s][j*s] == 2)
                    {
                        auto bv = cond.find_BV(Vec(j*h, i*h), Vec(j*h, (i-1)*h));
                        double alpha = (i*h - bv.p[1])/(bv.p[1] - (i-1)*h);
                        sum += -alpha*vec[to1d(i-1,j,n)];
                        tmp++;
                    }
                    vec[to1d(i,j,n)] = sum/double(tmp);
                }
            }
    }
    */
    if(_type == linear)
    {
        int s = originN / N;
        for(int i = 0 ; i <= N ; i ++)
            for(int j = 0 ; j <= N ; j ++)
            {
                if((!cond.isRegular) && (cond.loc[i*s][j*s] == 0))
                {
                    res(to1d(i,j,N)) = 0;
                    continue;
                }
                if(i%2 == 0)
                {
                    if(j%2 == 0) // x, y both even
                        res(to1d(i,j,N)) = vec(to1d(i/2,j/2,N/2));
                    else // x odd
                        res(to1d(i,j,N)) = (vec(to1d(i/2,(j+1)/2,N/2))+vec(to1d(i/2,(j-1)/2,N/2)))/2.0;
                }
                else
                {
                    if(j%2 == 0)
                        res(to1d(i,j,N)) = (vec(to1d((i+1)/2,j/2,N/2))+vec(to1d((i-1)/2,j/2,N/2)))/2.0;
                    else
                        res(to1d(i,j,N)) = (vec(to1d((i-1)/2,(j+1)/2,N/2))+vec(to1d((i-1)/2,(j-1)/2,N/2))+vec(to1d((i+1)/2,(j+1)/2,N/2))+vec(to1d((i+1)/2,(j-1)/2,N/2)))/4.0;
                }
            }
        return res;
    }
    else if(_type == quadratic) //TBD
    {
        ;
    }
    assert(0);
}

/*

template<>
VectorXd OP::restriction<1, full_weighting>(const VectorXd& _vec)
{
    const int N = _vec.size() - 1;
    assert(N % 2 == 0);
    VectorXd res(N/2 + 1);
    res(0) = _vec(0);
    res(N/2) = _vec(N);
    for(int j = 1 ; j < N/2 ; j ++)
        res[j] = (_vec[2*j-1]+2.0*_vec[2*j]+_vec[2*j+1])/4.0;
    return res;
}

template<>
VectorXd OP::restriction<2, full_weighting>(const VectorXd& _vec)
{
    int t = _vec.size();
    assert(isInteger(std::sqrt(t)));
    const int N = round(std::sqrt(t)) - 1;
    assert(N % 2 == 0);
    VectorXd res(SQUARE(N/2 + 1));
    

    for(int i = 0 ; i <= N/2 ; i ++)
        for(int j = 0 ; j <= N/2 ; j ++)
        {
            int s = 0;
            double sum = 0;
            if(i != 0)
            {
                sum += _vec(to1d(i*2,j*2,N))+_vec(to1d(i*2-1,j*2,N));
                s ++;
            }
            if(i != N/2)
            {
                sum += _vec(to1d(i*2,j*2,N))+_vec(to1d(i*2+1,j*2,N));
                s ++;
            }
            if(j != 0)
            {
                sum += _vec(to1d(i*2,j*2,N))+_vec(to1d(i*2,j*2-1,N));
                s ++;
            }
            if(j != N/2)
            {
                sum += _vec(to1d(i*2,j*2,N))+_vec(to1d(i*2,j*2+1,N));
                s ++;
            }
            res(to1d(i,j,N/2)) = sum/double(s*2);
        }
    return res;
}

template<>
VectorXd OP::restriction<1, injection>(const VectorXd& _vec)
{
    const int N = _vec.size() - 1;
    // std::cout << std::endl << N << std::endl;
    assert(N % 2 == 0);
    VectorXd res(N/2 + 1);
    for(int j = 0 ; j <= N/2 ; j ++)
        res[j] = _vec[2*j];
    return res;
}

template<>
VectorXd OP::restriction<2, injection>(const VectorXd& _vec)
{
    int t = _vec.size();
    assert(isInteger(std::sqrt(t)));
    const int N = round(std::sqrt(t)) - 1;
    assert(N % 2 == 0);
    VectorXd res(SQUARE(N/2 + 1));
    for(int i = 0 ; i <= N/2 ; i ++)
        for(int j = 0 ; j <= N/2 ; j ++)
            res[to1d(i,j,N/2)] = _vec[to1d(i*2,j*2,N)];
    return res;
}
// 0 1 2 3 4
// 0   1   2

template<>
VectorXd OP::interpolation<1, linear>(const VectorXd& _vec)
{
    const int N = (_vec.size()-1)*2;
    VectorXd res(N+1);
    for(int j = 0 ; j <= N ; j ++)
    {
        if(j%2 == 0)
            res[j] = _vec[j/2];
        else
            res[j] = (_vec[(j-1)/2]+_vec[(j+1)/2])/2.0;
    }
    return res;
}

template<>
VectorXd OP::interpolation<2, linear>(const VectorXd& _vec)
{
    int t = _vec.size();
    assert(isInteger(std::sqrt(t)));
    const int N = (round(std::sqrt(t)) - 1)*2;
    VectorXd res(SQUARE(N+1));
    for(int i = 0 ; i <= N ; i ++)
        for(int j = 0 ; j <= N ; j ++)
        {
            if(i%2 == 0)
            {
                if(j%2 == 0)
                    res(to1d(i,j,N)) = _vec(to1d(i/2,j/2,N/2));
                else
                    res(to1d(i,j,N)) = (_vec(to1d(i/2,(j+1)/2,N/2))+_vec(to1d(i/2,(j-1)/2,N/2)))/2.0;
            }
            else
            {
                if(j%2 == 0)
                    res(to1d(i,j,N)) = (_vec(to1d((i+1)/2,j/2,N/2))+_vec(to1d((i-1)/2,j/2,N/2)))/2.0;
                else
                    res(to1d(i,j,N)) = (_vec(to1d((i-1)/2,(j+1)/2,N/2))+_vec(to1d((i-1)/2,(j-1)/2,N/2))+_vec(to1d((i+1)/2,(j+1)/2,N/2))+_vec(to1d((i+1)/2,(j-1)/2,N/2)))/4.0;
            }
        }
    return res;
}

template<>
VectorXd OP::interpolation<1, quadratic>(const VectorXd& _vec) ////////TBD!!!
{
    const int N = (_vec.size()-1)*2;
    VectorXd res(N+1);
    for(int j = 0 ; j <= N ; j ++)
    {
        if(j%2 == 0)
            res[j] = _vec[j/2];
        else
            res[j] = (_vec[(j-1)/2]+_vec[(j+1)/2])/2.0;
    }
    return res;
}

template<>
VectorXd OP::interpolation<2, quadratic>(const VectorXd& _vec) ////////TBD!!!
{
    const int N = (_vec.size()-1)*2;
    VectorXd res(N+1);
    for(int j = 0 ; j <= N ; j ++)
    {
        if(j%2 == 0)
            res[j] = _vec[j/2];
        else
            res[j] = (_vec[(j-1)/2]+_vec[(j+1)/2])/2.0;
    }
    return res;
}


grid_operator get_restriction(int _dim, restriction_type _type)
{
    if(_dim == 1)
    {
        if(_type == full_weighting)
            return OP::restriction<1, full_weighting>;
        else if(_type == injection)
            return OP::restriction<1, injection>;
    }
    else if(_dim == 2)
    {
        if(_type == full_weighting)
            return OP::restriction<2, full_weighting>;
        else if(_type == injection)
            return OP::restriction<2, injection>;
    }
    assert(0);
}

grid_operator get_interpolation(int _dim, interpolation_type _type)
{
    if(_dim == 1)
    {
        if(_type == linear)
            return OP::interpolation<1, linear>;
        else if(_type == quadratic)
            return OP::interpolation<1, quadratic>;
    }
    else if(_dim == 2)
    {
        if(_type == linear)
            return OP::interpolation<2, linear>;
        else if(_type == quadratic)
            return OP::interpolation<2, quadratic>;
    }
    assert(0);
}
*/