/**
 * @file FD.h
 * @author Ouyang Shangke
 * @brief FD solver
 * @version 0.1
 * @date 2022-03-03
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#ifndef __FD_H__
#define __FD_H__

#include "Vec2d.h"
#include "conditions.h"
#include "global.h"
#include <Eigen/Sparse>

class FDsolver
{
public:
    Conditions cond; //< conditions of the Poisson equation
    Eigen::VectorXd res; //< store the result

public:
    /**
     * @brief Construct a new FDsolver object
     * 
     * @param _cond given conditions.
     */
    FDsolver(Conditions _cond);

    /**
     * @brief check if the condition is valid.
     * 
     * @return true for valid.
     * @return false for invalid.
     */
    bool check_cond();

    /**
     * @brief sovle the two-dimensional Poisson equation.
     * 
     * @return 0 for OK
     */
    int general_solver();

    /**
     * @brief Get the infty-, 1-, 2-norm errors.
     * 
     * @param u the exact solution
     * @param _type 0 for infty-norms; 1 for 1-norm and 2 for 2-norm;
     * @return the corresponding errors.
     */
    double get_err(double(*u)(double, double), int _type);

    /**
     * @brief Get the absoloute and relative errors on each discrete points,\n 
     * and output into a certain file (in the form of Ocatve code)
     * @param u the exact solution
     * @param _path path for output file
     * @param _index test index (default 0)
     * @return 0 for OK; -1 for err.
     */
    int get_err(double(*u)(double, double), std::string _path, int _index);

    /**
     * @brief print the solutions.
     * 
     * @param os outstream
     * @param _p given solver
     * @return
     */
    friend std::ostream &operator<<(std::ostream &os, const FDsolver &_p);
};

#endif