/**
 * @file FPN.h
 * @author Ouyang Shangke
 * @brief 
 * @version 0.1
 * @date 2021-12-06
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#ifndef __FPN_H__
#define __FPN_H__

#include "Config.h"

class FPN
{
private:
    /**
     * @brief Basic parameters for a FPN system
     * 
     */
    int beta, p, L, U;

public:

    /**
     * @brief Construct a new FPN object
     * 
     */
    FPN(){};    
    
    /**
     * @brief Construct a new FPN object
     * 
     * @param _beta 
     * @param _p 
     * @param _L 
     * @param _U 
     */
    FPN(int _beta, int _p, int _L, int _U): beta(_beta), p(_p), L(_L), U(_U){};

    /**
     * @brief get UFL of the FPN
     * 
     * @return UFL
     */
    double get_UFL();

    /**
     * @brief get OFL of the FPN
     * 
     * @return OFL
     */
    double get_OFL();

    /**
     * @brief Get all numbers in the normalized FPN system
     * 
     * @return a vector which stores demand numbers
     */
    vector_double get_all_norm();

    /**
     * @brief Get all the subnormal numbers in the FPN
     * 
     * @return a vector which stores demand numbers
     */
    vector_double get_all_subnorm();
};

#endif