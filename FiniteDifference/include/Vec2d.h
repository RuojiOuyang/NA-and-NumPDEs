/**
 * @file Vec2d.h
 * @author Ouyang Shangke
 * @brief Points in 2-dimention
 * @version 0.1
 * @date 2022-03-03
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#ifndef __VEC2D_H__
#define __VEC2D_H__

#include <ostream>
#include <map>
#include <math.h>
#include <initializer_list>

class Vec
{
protected:
    double coord[2];

public:
    Vec(){};

    Vec(double x, double y)
    {
        coord[0] = x;
        coord[1] = y;
    }
    // accessors
    double &operator[](int _d) { return coord[_d]; }
    const double &operator[](int _d) const { return coord[_d]; }

    #define ELMWISE_BINARY_OP(OpName, Op)                          \
    auto OpName(const Vec &rhs) const {                            \
        Vec res;                                                   \
        for (int i = 0; i < 2; i++) res[i] = coord[i] Op rhs[i];   \
        return res;                                                \
    }
    ELMWISE_BINARY_OP(operator+, +)
    ELMWISE_BINARY_OP(operator-, -)
    ELMWISE_BINARY_OP(operator*, *)
    ELMWISE_BINARY_OP(operator/, /)
    #undef ELMWISE_BINARY_OP

    #define RIGHT_BROADCAST(OpName, Op)                         \
    auto OpName(const double &rhs) const {                      \
        Vec res;                                                \
        for (int d = 0; d < 2; ++d) res[d] = coord[d] Op rhs;   \
        return res;                                             \
    }

    RIGHT_BROADCAST(operator+, +)
    RIGHT_BROADCAST(operator-, -)
    RIGHT_BROADCAST(operator*, *)
    RIGHT_BROADCAST(operator/, /)
    #undef RIGHT_BROADCAST

    int operator=(const Vec &p) {
        coord[0] = p[0];
        coord[1] = p[1];
        return 0;
    }

    bool operator==(const Vec &p) const{
        if(coord[0] == p[0] && coord[1] == p[1])
            return true;
        return false;
    }
    
    bool operator<(const Vec &p) const
    {
        if(coord[0] < p[0])
            return true;
        else if (coord[0] > p[0])
            return false;
        else
        {
            if(coord[1] < p[1])
                return true;
            else
                return false;
        }
    }

    double norm()
    {
        return std::sqrt(coord[0]*coord[0] + coord[1]*coord[1]);
    }

    friend std::ostream &operator<<(std::ostream &os, const Vec &p) {
        os << "(" << p[0] << " , " << p[1] << ")";
        return os;
    }

    friend double dist(Vec p1, Vec p2)
    {
        return std::sqrt((p1[0]-p2[0])*(p1[0]-p2[0]) + (p1[1]-p2[1])*(p1[1]-p2[1]));
    }
    
    friend double inner_product(Vec p1, Vec p2)
    {
        return p1[0]*p2[0] + p1[1]*p2[1];
    }
};

typedef std::map<Vec, double> FUNC_VALUE;

#endif