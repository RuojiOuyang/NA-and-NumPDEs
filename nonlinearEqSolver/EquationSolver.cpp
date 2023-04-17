#ifndef __EQUATIONSOLVER_CPP__
#define __EQUATIONSOLVER_CPP__

#include "EquationSolver.h"

double bisection::solve() const
{
    double left = this -> interval[0];
    double right = this -> interval[1];
    double mid = (left + right) / 2;
    double length = right - left;
    double flag = 0;

    for(int i = 1 ; i <= this -> M ; i ++)
    {
        if(std::abs(this -> f(mid)) < epsilon || length < delta)
        {
            return mid;
        }
        length = length / 2;
        if(this -> f(mid) * this -> f(left) <= 0)
            right = mid;
        else
            left = mid;
        mid = (right + left) / 2;
    }

    return D_ERR;
}

int bisection::set_interval(double _lowwer, double _upper)
{
    this->interval[0] = _lowwer;
    this->interval[1] = _upper;
    return 0;
}

double Newton::solve() const
{
    double flag = 0;
    double x = this->initial;
    for(int i = 1 ; i <= this->M ; i ++)
    {
        double t = this->f(x);
        if(std::abs(t) < epsilon)
        {
            return x;
        }
        x = x - t / this->df(x);
    }
    return D_ERR;
}

int Newton::set_initial(double _initial)
{
    this->initial = _initial;
    return 0;
}

double secant::solve() const
{
    double x_last = this->x0;
    double x_new = this->x1;
    double fx_new = this->f(x_new);
    double fx_last = this->f(x_last);
    
    for(int i = 1 ; i <= this->M ; i ++)
    {
        if(abs(fx_new) > abs(fx_last))
        {
            std::swap(x_new , x_last);
            std::swap(fx_new , fx_last);
        }
        double t = x_new - fx_new * ((x_new - x_last) / (fx_new - fx_last));
        x_last = x_new;
        fx_last = fx_new;
        x_new = t;
        fx_new = this->f(x_new);
        if(abs(x_new - x_last) < delta || abs(fx_new) < this->epsilon)
            return x_new;
    }

    return D_ERR;
}

int secant::set_initial(double _x0, double _x1)
{
    this->x0 = _x0;
    this->x1 = _x1;
    return 0;
}

#endif