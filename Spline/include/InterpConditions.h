/**
 * @file InterpConditions.h
 * @author OuyangShangke
 * @brief This file is for the class InterpConditions.
 * @version 0.1
 * @date 2021-11-23
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#ifndef __INTERPCONDITIONS_H__
#define __INTERPCONDITIONS_H__

#include "Config.h"
#include <iostream>

struct InterpConditions
{
    vector_double x; ///< interpolation sites
    vector2_double fx; ///< corresponding function values

    /**
     * @brief Default constructor.
     */
    InterpConditions(){};

    /**
     * @brief Construct a new Interp Conditions object.
     * 
     * @param list1 the initializer list of x
     * @param list2 the initializer list of fx
     */
    InterpConditions(std::initializer_list<double> list1, std::initializer_list<vector_double> list2)
    {
	    for(auto it = list1.begin(); it != list1.end(); it ++)
		    x.push_back(*it);
	    for(auto it = list2.begin(); it != list2.end(); it ++)
		    fx.push_back(*it);
    }
    /**
     * @brief Get the maximum order of the derivatives.
     * 
     * @return maximum order of the derivatives
     */
    int get_max_diff() const;
    /**
     * @brief Get the amount of all the function and its derivatives values.
     * 
     * @return total amount
     */
    int get_N() const;
    /**
     * @brief Check if the given point is one of the interpolation points.
     * 
     * @param _x given point
     * @return 
     *  @retval >=0 its position
     *  @retval -1 its not a interpolation point
     */
    int find_x(double _x) const;
    /**
     * @brief Input a InterpCondtions by a given structure.
     * 
     * @param input 
     * @param _interp given structure
     * @return
     */
    friend std::istream& operator>> (std::istream& input, InterpConditions& _interp);
};

#endif