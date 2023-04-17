/**
 * @file genJson.h
 * @author OuyangShangke
 * @brief help users to generate input json file
 * @version 0.1
 * @date 2022-03-25
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#ifndef __GENJSON_H__
#define __GENJSON_H__

#include "global.h"
#include "Vec2d.h"

template<int dim> class genJson;

/**
 * @brief generate json for 1-dimension problem
 * 
 * @tparam  
 */
template<>
class genJson<1>
{
private:
    int N;
    double1d f_val;
    double BV[2][2];

public:
    /**
     * @brief Construct a new gen Json object
     * 
     * @param _N the number of cells
     */
    genJson(int _N):N(_N){;};

    /**
     * @brief Set the -Laplace(u)
     * 
     * @param _f function
     * @return int 0 for OK
     */
    int set_fval(FUNC1d _f);

    /**
     * @brief set boundary value
     * 
     * @param _p point (0 or 1)
     * @param _type Dirichlet(0) or Neumann(1)
     * @param _val value
     * @return int 0 for OK
     */
    int set_BV(int _p, int _type, double _val);

    /**
     * @brief output to json
     * 
     * @param _path path for output json file
     * @return int 0 for OK
     */
    int write_json(std::string _path);
};

/**
 * @brief generate json for 2-dimension problem
 * 
 * @tparam  
 */
template<>
class genJson<2>
{
private:
    int N;
    double2d f_val;
    uint2d loc;
    bool isRegular;
    struct BoundaryValues //< store a single boundary value
    {
        Vec p; //< position of point
        Vec n_vec; //< direction of the derivative
        double val; //< value of function or its derivative
        BoundaryValues(Vec _p, Vec _n, double _val):p(_p),n_vec(_n),val(_val){;}; //< constructor
    };
    std::vector<BoundaryValues> BV;

public:
    /**
     * @brief Construct a new genJson object
     * 
     * @param _N the number of cells on each dimension
     * @param _isReg whether the domain is regular
     */
    genJson(int _N, bool _isReg):N(_N), isRegular(_isReg){;};

    /**
     * @brief set the json for sine domain. \n
     * This domain replaced the lower bound of (0,1)^2 with A*sin(pi*x), 
     * where A is the spectrum less than 1.
     * @param _spec spectrum 
     * @param _f -Laplace(u)
     * @param _u u(x,y)
     * @return int 0 for OK; -1 for error(invalid spectrum)
     */
    int set_sine_domain(double _spec, FUNC2d _f, FUNC2d _u);

    /**
     * @brief set the -Laplace(u). \n 
     * ONLY FOR REGULAR DOMAIN
     * @param _f -Laplace(u)
     * @return int 
     */
    int set_fval(FUNC2d _f);

    /**
     * @brief set the Dirichlet boundary value. \n
     * ONLY FOR REGULAR DOMAIN
     * @param _u u(x,y)
     * @return int 0 for OK
     */
    int set_BV0(FUNC2d _u);

    /**
     * @brief set the Neumann boundary value. \n
     * ONLY FOR REGULAR DOMAIN
     * @param _ux partial u/ partial x
     * @param _uy partial u/ partial y
     * @param _n direction vector
     * @return int 0 for OK
     */
    int set_BV1(FUNC2d _ux, FUNC2d _uy, Vec _n);

    /**
     * @brief set boundary value
     * 
     * @param _p points
     * @param _n direction vector. If it's Dirichlet funciton, then set it as (0,0)
     * @param _val value
     * @return int 0 for OK
     */
    int set_BV(Vec _p, Vec _n, double _val);

    /**
     * @brief output to json
     * 
     * @param _path path for output json file
     * @return int 0 for OK
     */
    int write_json(std::string _path);
};


// {
// private:    
//     double h; ///< grid size
//     Vec center; ///< center of circle
//     double radius; ///< raidus of circle
//     struct BoundaryValues //< Boundary Values
//     {
//         Vec p; ///< coordinate of point
//         Vec n_vec; ///< direction of derivatives
//         double val; ///< value
//         /**
//          * @brief Construct a new Boundary Values object
//          * 
//          * @param _p coordinate
//          * @param _n direction
//          * @param _val value
//          */
//         BoundaryValues(Vec _p, Vec _n, double _val):p(_p),n_vec(_n),val(_val){;};
//     };
//     double2d f_val; ///< store the RHS funciton values
//     std::vector<BoundaryValues> BV;

// public:
//     /**
//      * @brief Construct a new gen Json object
//      * 
//      * @param _h grid size 
//      */
//     genJson(double _h);

//     bool isInDomain(Vec _p);
//     bool isInSquare(Vec _p);
//     Vec toGrid(Vec _p);

//     /**
//      * @brief Set the circle domain 
//      * 
//      * @param center center of circle
//      * @param radius raidus of circle
//      * @return
//      *  @retval 0 for ok
//      *  @retval -1 for err (maintain connect, and cover at least 4 points)
//      */
//     int set_circ(Vec center, double radius);

//     /**
//      * @brief Set the RHS function values
//      * 
//      * @param f RHS function.
//      * @return 0 for OK
//      */
//     int set_fval(double (*f)(double, double));

//     /**
//      * @brief Set the fval object
//      * 
//      * @param _i location
//      * @param _j location
//      * @param _val value
//      */
//     void set_fval(int _i, int _j, double _val) { f_val[_i][_j] = _val; };

//     /**
//      * @brief set boundary values
//      * 
//      * @param _p coordinate
//      * @param _n direction
//      * @param _val values
//      */
//     void set_BV(Vec _p, Vec _n, double _val);
    
//     std::vector<Vec> discrete_circ();

//     /**
//      * @brief Set the Dirichlet conditions
//      * 
//      * @param u function
//      */
//     void set_BV0(double (*u)(double, double));

//     /**
//      * @brief Set the Dirichlet conditions for the boundary of square
//      * 
//      * @param u function
//      */
//     void set_square_BV0(double (*u)(double, double));

//     /**
//      * @brief Set the Dirichlet conditions for the boundary of circle
//      * 
//      * @param u function
//      */
//     void set_circle_BV0(double (*u)(double, double));

//     /**
//      * @brief Set the Neumann conditions
//      * 
//      * @param ux partial u/ partial x
//      * @param uy partial u/ partial y
//      * @param _n direction of derivative
//      */
//     void set_BV1(double (*ux)(double, double), double (*uy)(double, double), Vec _n);

//     /**
//      * @brief Set the Neumann conditions for the boundary of square
//      * 
//      * @param ux partial u/ partial x
//      * @param uy partial u/ partial y
//      * @param _n direction of derivative
//      */
//     void set_square_BV1(double (*ux)(double, double), double (*uy)(double, double), Vec _n);

//     /**
//      * @brief Set the Neumann conditions for the boundary of circle
//      * 
//      * @param ux partial u/ partial x
//      * @param uy partial u/ partial y
//      * @param _n direction of derivative
//      */
//     void set_circle_BV1(double (*ux)(double, double), double (*uy)(double, double), Vec _n);

//     /**
//      * @brief output to a json file
//      * 
//      * @param _path path of output json file
//      * @return 0 for OK
//      */
//     int write_json(std::string _path);
// };

#endif