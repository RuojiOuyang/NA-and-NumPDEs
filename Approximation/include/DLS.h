#ifndef __DLS_H__
#define __DLS_H__

#include "Config.h"

/**
 * @brief store points and function values.
 * 
 */
struct Conditions
{
    vector_double x; ///< points
    vector_double y; ///< function values
};

/**
 * @brief Discrete least square
 * 
 */
class DLS
{
private:
    double a0, a1, a2; ///< the quadratic polynomial p(x) = a2 x^2 + a1 x + a0

    static double f0(double _x){return 1;};
    static double f1(double _x){return _x;};
    static double f2(double _x){return _x*_x;};
    
public:
    /**
     * @brief Inner product of two functions.
     * 
     * @param _f1 func 1
     * @param _f2 func 2
     * @param _x points
     * @return result of inner product
     */
    friend double inner_product(double (*_f1)(double), double(*_f2)(double), vector_double _x);
    
    /**
     * @brief Inner product of a functions and a set.
     * 
     * @param _f1 func 1
     * @param _f set of function values
     * @param _x points
     * @return result of inner product
     */
    friend double inner_product(double (*_f1)(double), vector_double _f, vector_double _x);

    /**
     * @brief output a DLS object.
     * 
     * @param output 
     * @param _DLS DLS object
     * @return
     */
    friend std::ostream& operator<< (std::ostream& output, const DLS& _DLS);

    /**
     * @brief Solve DLS by normal equations.
     * 
     * @param _cond conditions
     * @return 0 for OK
     */
    int solve(const Conditions _cond);

    int solve_byQR(const Conditions _cond);
};

/**
 * @brief output a matrix which is stored in vector2_double
 * 
 * @param output 
 * @param _A matrix
 * @return
 */
std::ostream& operator<< (std::ostream& output, const vector2_double& _A);

#endif