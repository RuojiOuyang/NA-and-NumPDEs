/**
 * @file test_helper.h
 * @author OuyangShangke
 * @brief 
 * @version 0.1
 * @date 2022-04-16
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#ifndef __TEST_HELPER__
#define __TEST_HELPER__

#include <string>
#include <Eigen/Sparse>
#include "genJson.h"
#include "Multigrid.h"



using FUNC1d = double(*)(double);
using FUNC2d = double(*)(double, double);
using std::string, Eigen::VectorXd;

class helper
{
protected:
    string mark;
    std::vector<string> json_path;
    /**
     * type = 0 => Dirichlet
     * type = 1 => Mixed
     */
    int type;
    virtual int generate_json() = 0;
public:
    helper(int _t, const string &_mark):type(_t),mark(_mark){;};
    int output(string _title, string _filename, string _type, int _dim, VectorXd _x);
    virtual double get_abs_err(const VectorXd& _x, int _norm) = 0;
    virtual VectorXd get_abs_err(const VectorXd& _x) = 0;
    virtual int start_test() = 0;
};

class helper1d : public helper
{
private:
    FUNC1d u, f;
    double neumann;
    virtual int generate_json() override;
    // int test_all_type(Multigrid<1> &_mg);
public:
    helper1d(int _t, const string &_mark, FUNC1d _u, FUNC1d _f, double _n = 0)
        :helper(_t, _mark), u(_u), f(_f), neumann(_n){;};
    virtual double get_abs_err(const VectorXd& _x, int _norm) override; 
    virtual VectorXd get_abs_err(const VectorXd& _x) override;
    virtual int start_test() override;
};

class helper2d : public helper
{
private:
    FUNC2d u, f, ux, uy;
    Vec n_vec;
    double spectrum;
    virtual int generate_json() override;

public:
    helper2d(int _t, const string &_mark, FUNC2d _u, FUNC2d _f, double _spec = 0, FUNC2d _ux = NULL, FUNC2d _uy = NULL, Vec _n = Vec(0,0))
        :helper(_t, _mark), u(_u), f(_f), spectrum(_spec), ux(_ux), uy(_uy), n_vec(_n){;};
    virtual double get_abs_err(const VectorXd& _x, int _norm) override; 
    virtual VectorXd get_abs_err(const VectorXd& _x) override;
    virtual int start_test() override;
};

// class helper2d_irreg : public helper
// {
// private:
//     FUNC2d u, f, ux, uy;
//     virtual int generate_json() override;
//     // int test_all_type(Multigrid<1> &_mg);
// public:
//     helper1d(int _t, const string &_mark, FUNC1d _u, FUNC1d _f, double _n = 0)
//         :helper(_t, _mark), u(_u), f(_f), neumann(_n){;};
//     virtual double get_abs_err(const VectorXd& _x, int _norm) override; 
//     virtual VectorXd get_abs_err(const VectorXd& _x) override;
//     virtual int start_test() override;
// };


#endif