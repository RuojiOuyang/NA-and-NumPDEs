/**
 * @file global.h
 * @author Ouyang Shangke
 * @brief global func, constants, typedef, ...
 * @version 0.1
 * @date 2022-03-25
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#ifndef __PUBLIC_H__
#define __PUBLIC_H__

#include <math.h>
#include <vector>
#include <map>
#include <limits>
#include <string.h>
#include "Vec2d.h"

#define LEFT(point) (point-Vec(h,0))
#define RIGHT(point) (point+Vec(h,0))
#define UP(point) (point+Vec(0,h))
#define DOWN(point) (point-Vec(0,h))
#define SQUARE(x) ((x)*(x))

#define to1d(i,j) ((i)*(N+1)+(j))
#define to2d1(n) (((n)/(N+1)))
#define to2d2(n) (((n) - (N+1)*((n)/(N+1))))

// typedef
typedef std::vector<double> double1d;
typedef std::vector<std::vector<double> > double2d;
typedef std::string cstring;

// constant
const double EPS = 1e-8;
const double D_MAX = std::numeric_limits<double>::max(); ///< maximum of double

// func

bool isInteger(double t);

bool isInSquare(Vec _p);

double round(double t);

double max(double x, double y);

#endif