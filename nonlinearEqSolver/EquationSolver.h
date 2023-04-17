/**
 * @file EquationSolver.h
 * @author OuyangShangke
 * @brief 
 * @version 0.1
 * @date 2021-09-24
 * 
 * @copyright Copyright (c) 2021
 * 
 */
#ifndef __EQUATIONSOLVER_H__
#define __EQUATIONSOLVER_H__

#include <iostream>
#include <math.h>
#include <limits>
#include <iomanip>
#include <vector>
#include <chrono>

#define D_MAX (std::numeric_limits<double>::max())
#define D_MIN (std::numeric_limits<double>::min())
#define D_ERR (std::numeric_limits<double>::max()-1)

#define default_epsilon 1e-8
#define default_delta 1e-5
#define default_M 1e5

/**
 * @brief 
 * 
 */
class EquationSolver
{
protected:
    double (*f)(double); ///< function
    double epsilon;      ///< accuracy of solution
    double M;            ///< limit of iteration

    /**
     * @brief Pure virtual function.
     * 
     * @return double 
     */
    virtual double solve() const = 0;

    /**
     * @brief Construct a new Equation Solver object
     * 
     * @param _f given function
     * @param _M limit of iteration
     * @param _epsilon accuracy of solution
     */
    EquationSolver(double (*_f)(double), double _M = default_M, double _epsilon = default_epsilon) 
    : f(_f), M(_M), epsilon(_epsilon){};
};

/**
 * @brief Bisection method class.
 * Derived from EquationSolver class
 */
class bisection : public EquationSolver
{
private:
    double interval[2];
    double delta; ///< limits of the interval length

public:
    /**
     * @brief Construct a new bisection object
     * 
     * @param _f given function
     * @param _lower lowwer bound of the interval
     * @param _upper upper bound of the interval
     * @param _M limit of iteration
     * @param _epsilon accuracy of the solution
     * @param _delta limit of the interval length
     */
    bisection(double (*_f)(double), double _lower, double _upper, double _M = default_M, double _epsilon = default_epsilon, 
        double _delta = default_delta)
        :EquationSolver(_f, _M, _epsilon), delta(_delta), interval{_lower , _upper}{};

    /**
     * @brief Reset the searching interval.
     * 
     * @param _lowwer lowwer bound of the interval
     * @param _upper upper bound of the interval
     * @return int 0 for OK
     */
    int set_interval(double _lowwer, double _upper);

    /**
     * @brief Bisection method.
     * 
     * @return double solution
     */
    double solve() const;
};

/**
 * @brief Newton's method class.
 * Derived from EquationSolver class
 */
class Newton : public EquationSolver
{
private:
    double (*df)(double); ///< derivative of f(x)
    double initial; ///< initial point of Newton's Method

public:
    /**
     * @brief Construct a new Newton object
     * 
     * @param _f given function
     * @param _df the derivative of the given function
     * @param _initial initial point of Newton iteration
     * @param _M limit of iteration
     * @param _epsilon accuracy of the solution
     */
    Newton(double (*_f)(double), double (*_df)(double), double _initial, double _M = default_M, double _epsilon = default_epsilon)
        :EquationSolver(_f, _M, _epsilon), df(_df), initial(_initial){};

    /**
     * @brief Reset initial point of Newton's method.
     * 
     * @param _initial initial point
     * @return int 0 for OK
     */
    int set_initial(double _initial);

    /**
     * @brief Newton's method.
     * 
     * @return double solution
     */
    double solve() const;
};

/**
 * @brief Secant method class.
 * Derived from EquationSolver class
 */
class secant : public EquationSolver
{
private:
    double x0; ///< initial point x0
    double x1; ///< initial point x1
    double delta; ///< limit of the distance of 2 points from adjacent iteration

public:
    /**
     * @brief Construct a new secant object
     * 
     * @param _f given function
     * @param _x0 initial point x0
     * @param _x1 initial point x1
     * @param _M limit of iteration
     * @param _epsilon accuracy of the solution
     * @param _delta limit of the distance of 2 points from adjacent iteration
     */
    secant(double (*_f)(double), double _x0, double _x1, double _M = default_M, double _epsilon = default_epsilon, double _delta = default_delta)
        :EquationSolver(_f, _M, _epsilon), x0(_x0), x1(_x1), delta(_delta){};

    /**
     * @brief Reset initial points.
     * 
     * @param _x0 initial point x0
     * @param _x1 initial point x1
     * @return int 0 for OK 
     */
    int set_initial(double _x0, double _x1);

    /**
     * @brief secant method.
     * 
     * @return double solution
     */
    double solve() const;
};

#endif