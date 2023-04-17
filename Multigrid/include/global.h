/**
 * @file global.h
 * @author OuyangShangke
 * @brief global func, constants, typedef, ...
 * @version 0.1
 * @date 2022-03-25
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#ifndef __GLOBAL_H__
#define __GLOBAL_H__

#include <iostream>
#include <fstream>
#include <sstream>
#include <json/json.h>
#include <assert.h>
#include <math.h>
#include <vector>
#include <map>
#include <limits>
#include <string.h>
#include <Eigen/Sparse>

#define LEFT(point) (point-Vec(h,0))
#define RIGHT(point) (point+Vec(h,0))
#define UP(point) (point+Vec(0,h))
#define DOWN(point) (point-Vec(0,h))
#define SQUARE(x) ((x)*(x))

// #define to1d(i,j) ((i)*(N+1)+(j))
#define to2d1(n) (((n)/(N+1)))
#define to2d2(n) (((n) - (N+1)*((n)/(N+1))))

// typedef
using double1d = std::vector<double>;
using double2d = std::vector<double1d>;
using int1d = std::vector<int>;
using int2d = std::vector<std::vector<int> >;
using uint1d = std::vector<unsigned int>;
using uint2d = std::vector<std::vector<unsigned int> >;
using string = std::string;
using FUNC1d = double (*)(double);
using FUNC2d = double (*)(double, double);

// in Eigen
using VectorXd = Eigen::VectorXd;
using SparseMatrix = Eigen::SparseMatrix<double>;
using grid_operator = VectorXd(*)(const VectorXd&);
// using double1d =  std::vector<double>;
// using double2d =  std::vector<std::vector<double> >;
// using int1d = std::vector<int>;
// using int2d = std::vector<std::vector<int> >;
// using uint1d = std::vector<unsigned int>;
// using uint2d = std::vector<std::vector<unsigned int> >;

// typedef std::string cstring;


// constant
const double EPS = 1e-8;
const double D_MAX = std::numeric_limits<double>::max(); ///< maximum of double

// func

bool isInteger(double t);

bool isEqual(double a, double b);

double round(double t);

double max(double x, double y);

int to1d(int i , int  j , int N);

#endif