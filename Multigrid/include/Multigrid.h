/**
 * @file Multigrid.h
 * @author OuyangShangke
 * @brief 
 * @version 0.1
 * @date 2022-04-04
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#ifndef __MULTIGRID_H__
#define __MULTIGRID_H__

#include "global.h"
#include "operator.h"
#include "Generator.h"

template<int dim>
class Multigrid
{
public:
    int MAX_iteration; 
    double epsilon;
    restriction_type res_type;
    interpolation_type inter_type;
    // grid_operator restriction;
    // grid_operator interpolation;
    std::vector<SparseMatrix> Matrix;
    VectorXd f;
    OP<dim>* op;
    // uint2d loc;
    double1d VC_recorder;

public:

    // Multigrid(int _M, double _eps, restriction_type _res, interpolation_type _inter, const std::vector<SparseMatrix> &_Matrix, const VectorXd &_f)
    //     :MAX_iteration(_M), epsilon(_eps), myRestriction(_res), myInterpolation(_inter), Matrix(_Matrix), f(_f){;};

    Multigrid(int _M, double _eps, restriction_type _res, interpolation_type _inter, const string &_path);

    // Multigrid(int _M, double _eps, restriction_type _res, interpolation_type _inter, Generator<dim> &_Generator);

    void set_MAX_iteration(int _M){MAX_iteration = _M;};

    void set_epsilon(double _eps){epsilon = _eps;};

    void set_restriction(restriction_type _res){res_type = _res;};

    void set_interpolation(interpolation_type _inter){inter_type = _inter;};

    VectorXd FMG_solver(int _nu1, int _nu2);

    VectorXd FMG_solver(VectorXd _f, int _nu1, int _nu2, int _depth = 0);

    VectorXd VC_solver(int _nu1, int _nu2, bool _recorder_on = false, VectorXd _initial = VectorXd());

    VectorXd VC_solver(VectorXd _initial, VectorXd _f, int _nu1, int _nu2, int _depth = 0);

    VectorXd LU_solver();

    double get_abs_err(const VectorXd& _x, int _t = 0);

    double get_rel_err(const VectorXd& _x, int _t = 0);

    double1d get_recorder(){return VC_recorder;};
};

VectorXd relax(SparseMatrix _A, VectorXd _b, VectorXd _initial, int _nu, double _omega = 2.0/3.0);

#endif