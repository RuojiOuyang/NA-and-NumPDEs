/**
 * @file Multigrid.cpp
 * @author OuyangShangke
 * @brief 
 * @version 0.1
 * @date 2022-04-04
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include "Multigrid.h"

template<>
Multigrid<1>::Multigrid(int _M, double _eps, restriction_type _res, interpolation_type _inter, const string &_path)
    :MAX_iteration(_M), epsilon(_eps), inter_type(_inter), res_type(_res)
{
    // restriction = get_restriction(1, _res);
    // interpolation = get_interpolation(1, _inter);
    Conditions<1> cond(_path);
    Generator<1> gen(cond);
    gen.generate(Matrix, f);
    op = new OP<1>;
}

template<>
Multigrid<2>::Multigrid(int _M, double _eps, restriction_type _res, interpolation_type _inter, const string &_path)
    :MAX_iteration(_M), epsilon(_eps), inter_type(_inter), res_type(_res)
{
    // restriction = get_restriction(2, _res);
    // interpolation = get_interpolation(2, _inter);
    Conditions<2> cond(_path);
    Generator<2> gen(cond);
    gen.generate(Matrix, f);
    
    op = new OP<2>(cond);
    // if(!cond.isRegular)
    // {
    //     // std::cout << "sss" << std::endl;
        
    // }
    // else
    // {
    //     // std::cout << "ggg" << std::endl;
    //     op = new OP<2>(cond.N);
    // }
}

// template<int dim>
// Multigrid<dim>::Multigrid(int _M, double _eps, restriction_type _res, interpolation_type _inter, Generator<dim> &_Generator)
//     :MAX_iteration(_M), epsilon(_eps)
// {
//     restriction = get_restriction(dim, _res);
//     interpolation = get_interpolation(dim, _inter);
//     _Generator.generate(Matrix, f);
// }

template<int dim>
VectorXd Multigrid<dim>::VC_solver(int _nu1, int _nu2, bool _recorder_on, VectorXd _initial)
{
    int N = f.size();
    double1d recorder;
    if(_initial.size() == 0)
        _initial = VectorXd::Zero(N);
    assert(_initial.size() == N);
    auto &v = _initial;
    // std::cout << "_________________" << std::endl;
    for(int i = 1 ; i <= MAX_iteration ; i ++)
    {
        double tmp = get_rel_err(v);
        recorder.push_back(tmp);
        if(tmp <= epsilon)
            break;
        // std::cout << tmp << std::endl;
        v = VC_solver(v, f, _nu1, _nu2);
    }
    if(_recorder_on)
        VC_recorder = recorder;
    return v;
}

template<int dim>
VectorXd Multigrid<dim>::VC_solver(VectorXd _initial, VectorXd _f, int _nu1, int _nu2, int _depth)
{
    // relax
    // std::cout << "depth1: " << _depth << std::endl;
    if(_initial.size() == 0)
        _initial = VectorXd::Zero(_f.size());

    const int N = _initial.size();
    VectorXd v = relax(Matrix[_depth], _f, _initial, _nu1);
    
    if(_depth == Matrix.size()-1) // coarsest
    {
        v = relax(Matrix[_depth], _f, v, _nu2);
        return v;
    }
    // std::cout << "depth2: " << _depth << std::endl;
    VectorXd v2 = VC_solver(VectorXd(), op->restriction(res_type, _f-Matrix[_depth]*v), _nu1, _nu2, _depth+1);
    v = v + op->interpolation(inter_type, v2);
    v = relax(Matrix[_depth], _f, v, _nu2);
    return v;
}

template<int dim>
VectorXd Multigrid<dim>::FMG_solver(int _nu1, int _nu2)
{
    return FMG_solver(f, _nu1, _nu2);
}

template<int dim>
VectorXd Multigrid<dim>::FMG_solver(VectorXd _f, int _nu1, int _nu2, int _depth)
{
    const int N = _f.size();
    VectorXd v;
    if(_depth == Matrix.size()-1)
    {
        v = Eigen::VectorXd::Zero(N);
    }
    else
    {
        VectorXd v2 = FMG_solver(op->restriction(res_type,_f), _nu1, _nu2, _depth+1);
        v = op->interpolation(inter_type,v2);
    }
    v = VC_solver(v, _f, _nu1, _nu2, _depth);
    return v;
}

template<int dim>
VectorXd Multigrid<dim>::LU_solver()
{
    SparseMatrix A = Matrix[0];
    A.makeCompressed();
    Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > solver;
    solver.analyzePattern(A);
    solver.factorize(A);
    return solver.solve(f);
}

template<int dim>
double Multigrid<dim>::get_abs_err(const VectorXd& _x, int _t)
{
    assert(_x.size() == f.size());
    VectorXd e = Matrix[0]*_x - f;
    if(_t == 0)
        return e.lpNorm<Eigen::Infinity>();
    if(_t == 1)
        return e.lpNorm<1>();
    if(_t == 2)
        return e.lpNorm<2>();
    return D_MAX;
}

template<int dim>
double Multigrid<dim>::get_rel_err(const VectorXd& _x, int _t)
{
    assert(f.lpNorm<1>() != 0 && _x.size() == f.size());
    VectorXd e = Matrix[0]*_x - f;
    if(_t == 0)
        return e.lpNorm<Eigen::Infinity>()/f.lpNorm<Eigen::Infinity>();
    if(_t == 1)
        return e.lpNorm<1>()/f.lpNorm<1>();
    if(_t == 2)
        return e.lpNorm<2>()/f.lpNorm<2>();
    return D_MAX;
}

VectorXd relax(SparseMatrix _A, VectorXd _b, VectorXd _initial, int _nu, double _omega)
{
    long N = _b.size();
    // std::cout << "relax: " << _A.size() << " " << _initial.size() << " " << _b.size() << " " << _initial.size() << std::endl;
    assert(_A.size()==N*N && _initial.size()==N);
    SparseMatrix diag(N, N);
    SparseMatrix T(N, N);
    VectorXd c(N);
    VectorXd &v = _initial;
    for(int i = 0 ; i < N ; i ++)
        diag.coeffRef(i,i) = 1/_A.coeffRef(i,i);
    _A = -_omega * diag * _A;
    for(int i = 0 ; i < N ; i ++)
        _A.coeffRef(i,i) += 1;
    c = _omega*diag*_b;
    for(int i = 1 ; i <= _nu ; i ++)
        v = _A*v + c;
    // std::cout << "after: " << v.size() << std::endl;
    return v;
}

template class Multigrid<1>;
template class Multigrid<2>;
