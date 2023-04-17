/**
 * @file conditions.h
 * @author OuyangShangke
 * @brief store the conditions of Poisson equation
 * @version 0.1
 * @date 2022-03-25
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#ifndef __CONDITIONS_H__
#define __CONDITIONS_H__

#include "global.h"
#include "Vec2d.h"

template<int dim> class Conditions; 

template<>
class Conditions<1>
{
public:
    int N;
    double1d f_val;
    double BV[2][2];

    Conditions(std::string _path);
};

template<>
class Conditions<2>
{
public:
    int N;
    double2d f_val;
    uint2d loc; // 0 - out domain, 1 - on boundary, 2 - in domain
    bool isRegular;
    struct BoundaryValues //< store a single boundary value
    {
        Vec p; //< position of point
        Vec n_vec; //< direction of the derivative
        double val; //< value of function or its derivative
        BoundaryValues(Vec _p, Vec _n, double _val):p(_p),n_vec(_n),val(_val){;}; //< constructor
    };
    std::vector<BoundaryValues> BV;

    Conditions(std::string _path);

    /**
     * @brief get the boundary information of a given point.
     * 
     * @param _p the given point
     * @return struct BoundaryValues& 
     */
    const struct BoundaryValues &find_BV(Vec _p);

    /**
     * @brief get the boundary points between p1 and p2. \n 
     * if there're several boundary points between p1 and p2, then we choose the one closed to p1 most.
     * @param _p1 given point 1
     * @param _p2 given point 2
     * @return struct BoundaryValues&
     */
    const struct BoundaryValues &find_BV(Vec _p1, Vec _p2);
};
#endif