/**
 * @file operator.h
 * @author OuyangShangke
 * @brief 
 * @version 0.1
 * @date 2022-04-14
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#ifndef __OPERATOR_H__
#define __OPERATOR_H__

#include "global.h"
#include "conditions.h"

enum restriction_type{full_weighting, injection};
enum interpolation_type{linear, quadratic};

template<int dim> class OP;

template<>
class OP<1>
{
public:
    VectorXd restriction(restriction_type _type, const VectorXd& _vec);
    VectorXd interpolation(interpolation_type _type, const VectorXd& _vec);
};

template<>
class OP<2>
{
private:
    // uint2d loc;
    Conditions<2> cond;
    bool isRegular;
public:
    VectorXd restriction(restriction_type _type, const VectorXd& _vec);
    VectorXd interpolation(interpolation_type _type, const VectorXd& _vec);
    OP(Conditions<2> &_cond):cond(_cond){;};
    // OP(int _N, Conditions<2> &_cond):N_(_N), cond(_cond), isRegular(false){;};
};

/**
 * Restriction
 */

// grid_operator get_restriction(int _dim, restriction_type _type);

// grid_operator get_interpolation(int _dim, interpolation_type _type);

// VectorXd res_d1_FW(const VectorXd& _vec);

// VectorXd res_d1_injection(const VectorXd& _vec);

// VectorXd inter_d1_linear(const VectorXd& _vec);

// VectorXd inter_d1_q(const VectorXd& _vec);
// namespace OP
// {
//     template<int dim, restriction_type type>
//     VectorXd restriction(const VectorXd& _vec);

//     template<int dim, interpolation_type type>
//     VectorXd interpolation(const VectorXd& _vec);
// }


#endif