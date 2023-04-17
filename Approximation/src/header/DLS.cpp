#include "DLS.h"
#include "myLapack.h"

double inner_product(double (*_f1)(double), double(*_f2)(double), vector_double _x)
{
    double sum = 0;
    for(auto it = _x.begin() ; it < _x.end() ; it ++)
        sum += _f1(*it) * _f2(*it);
    return sum;
}

double inner_product(double (*_f1)(double), vector_double _f, vector_double _x)
{
    double sum = 0;
    auto it_f = _f.begin();
    for(auto it = _x.begin(); it < _x.end() ; it ++, it_f++)
        sum += _f1(*it) * (*it_f);
    return sum;
}

int DLS::solve(const Conditions _cond)
{
    vector2_double A;
    auto &x = _cond.x;

    vector_double b = {inner_product(f0, _cond.y,x), inner_product(f1, _cond.y,x), inner_product(f2, _cond.y,x)};
    A.push_back(vector_double{inner_product(f0,f0,x),inner_product(f0,f1,x),inner_product(f0,f2,x)});
    A.push_back(vector_double{inner_product(f1,f0,x),inner_product(f1,f1,x),inner_product(f1,f2,x)});
    A.push_back(vector_double{inner_product(f2,f0,x),inner_product(f2,f1,x),inner_product(f2,f2,x)});
    solve_general(A, b);
    this->a0 = b[0];
    this->a1 = b[1];
    this->a2 = b[2];
    return 0;
}
int DLS::solve_byQR(const Conditions _cond)
{
    vector2_double A, Q, R;
    vector2_double Gram;
    vector_double b, tmp = {0, 0, 0};
    auto &x = _cond.x;
    auto &y = _cond.y;

    /// original Gram matrix
    Gram.push_back(vector_double{inner_product(f0,f0,x),inner_product(f0,f1,x),inner_product(f0,f2,x)});
    Gram.push_back(vector_double{inner_product(f1,f0,x),inner_product(f1,f1,x),inner_product(f1,f2,x)});
    Gram.push_back(vector_double{inner_product(f2,f0,x),inner_product(f2,f1,x),inner_product(f2,f2,x)});
    std::cout << "Gram:" << std::endl;
    std::cout << Gram << std::endl;

    for(auto it = y.begin() ; it < y.end() ; it++)
        b.push_back(*it);
    
    for(auto it = x.begin() ; it < x.end() ; it++)
        A.push_back(vector_double{f0(*it), f1(*it), f2(*it)});

    QR(A, Q, R);
    
    for(int i = 0 ; i < 3; i ++)
        for(int j = 0 ; j < y.size() ; j ++)
            tmp[i] += Q[j][i] * b[j];
    b.resize(3);

    for(int i = 0 ; i < 3 ; i ++)
        b[i] = tmp[i];
    for(int j = 3 ; j < y.size() ; j ++)
        R.pop_back();
    std::cout << "R:" << std::endl;
    std::cout << R << std::endl;
    
    solve_general(R, b);
    
    this->a0 = b[0];
    this->a1 = b[1];
    this->a2 = b[2];
    return 0;

}

std::ostream& operator<< (std::ostream& output , const DLS& _DLS)
{
    output << _DLS.a0 << " " << _DLS.a1 << " " << _DLS.a2 << std::endl;
    return output;
}

std::ostream& operator<< (std::ostream& output, const vector2_double& _A)
{
    for(int i = 0 ; i < _A.size() ; i ++)
    {
        for(int j = 0 ; j < _A[i].size() ; j ++)
            std::cout << _A[i][j] << " ";
        std::cout << std::endl;
    }
    return output;
}