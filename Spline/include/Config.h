/**
 * @file Config.h
 * @author OuyangShangke
 * @brief This file is for declarations of constants, macro definitions, and some public headers.
 * @version 0.1
 * @date 2021-11-23
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#ifndef __CONFIG_H__
#define __CONFIG_H__

/**
 * @brief Header
 * 
 */
#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <limits>
#include <math.h>
#include <initializer_list>
#include "Vec.h"

// CONSTANT
const double SCANNING_STEP = 1e-3; ///< the default step of scanning
const double EPSILON = 1e-8; ///< rough estimation of machine precision

// int
typedef unsigned int u_int; ///< unsigned int

// vector
typedef std::vector<double> vector_double; ///< vector for double
typedef std::vector<std::vector<double>> vector2_double; ///< vector for vector<double>
typedef std::vector<int> vector_int; ///< vector for int
typedef std::vector<std::string> vector_string; ///< vector for string

// map
typedef std::map<int, double> M_INT_DOUBLE; ///< map int -> double
typedef std::map<double, double> M_DOUBLE_DOUBLE; ///< map double -> double
typedef std::pair<int, double> P_INT_DOUBLE; ///< pair int -> double
typedef std::pair<double, double> P_DOUBLE_DOUBLE; ///< pair double -> double

// string
typedef std::string STRING; ///< std::string
const std::string NULL_str = STRING(""); ///< empty string

/// Vec
#define Vec1 Vec<double,1> ///< Vec for dimention 1
#define Vec2 Vec<double,2> ///< Vec for dimention 2
#define vector_Vec1 std::vector<Vec<double,1>> ///< vector for d-1 Vec
#define vector_Vec2 std::vector<Vec<double,2>> ///< vector for d-2 Vec

/// limits
#define D_MAX (std::numeric_limits<double>::max()) ///< maximum of double
#define D_ERR (std::numeric_limits<double>::max()-1) ///< symbol of err

#endif