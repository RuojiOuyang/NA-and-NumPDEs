/**
 * @file conditions.h
 * @author Ouyang Shangke
 * @brief store the conditions of Poisson equation
 * @version 0.1
 * @date 2022-03-25
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#ifndef __CONDITIONS_H__
#define __CONDITIONS_H__

#include <jsoncpp/json/json.h>
#include "global.h"
#include "Vec2d.h"

class Conditions
{
public:    
    double h; //< grid size
    Vec center; //< center of circle
    double radius; //< radius of circle
    struct BoundaryValues //< store a single boundary value
    {
        Vec p; //< position of point
        Vec n_vec; //< direction of the derivative
        double val; //< value of function or its derivative
        BoundaryValues(Vec _p, Vec _n, double _val):p(_p),n_vec(_n),val(_val){;}; //< constructor
    };
    double2d f_val; //< store the value of -Laplace(u)
    std::vector<BoundaryValues> BV; //< store the value of boundary values

    /**
     * @brief Construct a new Conditions object by a json.
     * 
     * @param _path path of the json
     */
    Conditions(std::string _path);

    double get_h() const{ return h; };

    int get_N() const{ return round(1/h); };
    
    struct BoundaryValues &find_BV(Vec _p);

    Vec toGrid(Vec _p);

    int find_i(int i, double x);

    int find_j(int j, double y);
};

/**
 条件类型：
 1、f，都在格点上
 2、边界条件
    1、常规情况在各点上（regular）
    2、非常规可能在特殊点上
 */


#endif