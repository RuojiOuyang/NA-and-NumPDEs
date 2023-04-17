/**
 * @file genJson.h
 * @author Ouyang Shangke
 * @brief help users to generate input json file
 * @version 0.1
 * @date 2022-03-25
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#ifndef __GENJSON_H__
#define __GENJSON_H__

#include <jsoncpp/json/json.h>
#include "global.h"
#include "vector"

class genJson
{
private:    
    double h; ///< grid size
    Vec center; ///< center of circle
    double radius; ///< raidus of circle
    struct BoundaryValues //< Boundary Values
    {
        Vec p; ///< coordinate of point
        Vec n_vec; ///< direction of derivatives
        double val; ///< value
        /**
         * @brief Construct a new Boundary Values object
         * 
         * @param _p coordinate
         * @param _n direction
         * @param _val value
         */
        BoundaryValues(Vec _p, Vec _n, double _val):p(_p),n_vec(_n),val(_val){;};
    };
    double2d f_val; ///< store the RHS funciton values
    std::vector<BoundaryValues> BV;

public:
    /**
     * @brief Construct a new gen Json object
     * 
     * @param _h grid size 
     */
    genJson(double _h);

    bool isInDomain(Vec _p);
    bool isInSquare(Vec _p);
    Vec toGrid(Vec _p);

    /**
     * @brief Set the circle domain 
     * 
     * @param center center of circle
     * @param radius raidus of circle
     * @return
     *  @retval 0 for ok
     *  @retval -1 for err (maintain connect, and cover at least 4 points)
     */
    int set_circ(Vec center, double radius);

    /**
     * @brief Set the RHS function values
     * 
     * @param f RHS function.
     * @return 0 for OK
     */
    int set_fval(double (*f)(double, double));

    /**
     * @brief Set the fval object
     * 
     * @param _i location
     * @param _j location
     * @param _val value
     */
    void set_fval(int _i, int _j, double _val) { f_val[_i][_j] = _val; };

    /**
     * @brief set boundary values
     * 
     * @param _p coordinate
     * @param _n direction
     * @param _val values
     */
    void set_BV(Vec _p, Vec _n, double _val);
    
    std::vector<Vec> discrete_circ();

    /**
     * @brief Set the Dirichlet conditions
     * 
     * @param u function
     */
    void set_BV0(double (*u)(double, double));

    /**
     * @brief Set the Dirichlet conditions for the boundary of square
     * 
     * @param u function
     */
    void set_square_BV0(double (*u)(double, double));

    /**
     * @brief Set the Dirichlet conditions for the boundary of circle
     * 
     * @param u function
     */
    void set_circle_BV0(double (*u)(double, double));

    /**
     * @brief Set the Neumann conditions
     * 
     * @param ux partial u/ partial x
     * @param uy partial u/ partial y
     * @param _n direction of derivative
     */
    void set_BV1(double (*ux)(double, double), double (*uy)(double, double), Vec _n);

    /**
     * @brief Set the Neumann conditions for the boundary of square
     * 
     * @param ux partial u/ partial x
     * @param uy partial u/ partial y
     * @param _n direction of derivative
     */
    void set_square_BV1(double (*ux)(double, double), double (*uy)(double, double), Vec _n);

    /**
     * @brief Set the Neumann conditions for the boundary of circle
     * 
     * @param ux partial u/ partial x
     * @param uy partial u/ partial y
     * @param _n direction of derivative
     */
    void set_circle_BV1(double (*ux)(double, double), double (*uy)(double, double), Vec _n);

    /**
     * @brief output to a json file
     * 
     * @param _path path of output json file
     * @return 0 for OK
     */
    int write_json(std::string _path);
};

#endif