#include "Generator.h"

using T = Eigen::Triplet<double>;
using vecT = std::vector<T>;

int Generator<1>::generate(std::vector<SparseMatrix>& _Matrix, VectorXd &_f)
{
    int N = cond.N;
    _f.resize(N+1);
    
    for(int i = 0 ; i <= N ; i ++) // generate vector f
    {
        if(i == 0)
            _f(i) = cond.BV[0][0];
        else if(i == N)
            _f(i) = cond.BV[1][0];
        else
            _f(i) = cond.f_val[i];
    }
    while(N >= coarsest_N)
    {
        SparseMatrix A(N+1, N+1);
        vecT tripletlist;
        double h = 1/double(N);
        for(int i = 0 ; i <= N ; i ++)
        {
            if(i == 0)
            {
                if(cond.BV[0][1] == 0)
                {
                    tripletlist.push_back(T(0, 0, 1));
                }
                else // Neumann
                {
                    tripletlist.push_back(T(0, 0, -3/2.0/h));
                    tripletlist.push_back(T(0, 1, 2/h));
                    tripletlist.push_back(T(0, 2, -1/2.0/h));
                }
            }
            else if(i == N)
            {
                if(cond.BV[1][1] == 0)
                    tripletlist.push_back(T(N, N, 1));
                else // Neumann
                {
                    tripletlist.push_back(T(N, N, 3/2.0/h));
                    tripletlist.push_back(T(N, N-1, -2/h));
                    tripletlist.push_back(T(N, N-2, 1/2.0/h));
                }
            }
            else
            {
                tripletlist.push_back(T(i,i-1, -1/SQUARE(h)));
                tripletlist.push_back(T(i,i,2/SQUARE(h)));
                tripletlist.push_back(T(i,i+1, -1/SQUARE(h)));
            }
        }
        A.setFromTriplets(tripletlist.begin(), tripletlist.end());
        _Matrix.push_back(A);
        N /= 2;
    }
    return 0;
}

int Generator<2>::generate(std::vector<SparseMatrix>& _Matrix, VectorXd &_f)
{
    if(cond.isRegular)
        gen_regular(_Matrix, _f);
    else
        gen_irregular(_Matrix, _f);
    return 0;
}

int Generator<2>::gen_regular(std::vector<SparseMatrix>& _Matrix, VectorXd &_f)
{
    int originN = cond.N;
    int N = originN;
    _f.resize(SQUARE(N+1));
    
    while(N >= coarsest_N)
    {
        SparseMatrix A(SQUARE(N+1), SQUARE(N+1));
        vecT tripletlist;
        double h = 1/double(N);

        for(int i = 0 ; i <= N ; i ++)
            for(int j = 0 ; j <= N ; j ++)
            {
                if(i == 0 || j == 0 || i == N || j == N)
                {
                    auto bv = cond.find_BV(Vec(j*h, i*h));
                    if(N == originN)
                        _f(to1d(i,j,N)) = bv.val;
                    if(bv.n_vec == Vec(0,0)) // Dirichlet
                    {
                        tripletlist.push_back(T(to1d(i,j,N), to1d(i,j,N), 1));
                        continue;
                    }
                    // Neumann
                    double center_coeff = 0;
                    double x = bv.n_vec[0], y = bv.n_vec[1];
                    if(i == 0) // y direction
                    {
                        center_coeff += -3/2.0/h*y;
                        tripletlist.push_back(T(to1d(i,j,N), to1d(i+1,j,N), 2.0/h*y));
                        tripletlist.push_back(T(to1d(i,j,N), to1d(i+2,j,N), -1/2.0/h*y));
                    }
                    else if(i == N)
                    {
                        center_coeff += 3/2.0/h*y;
                        tripletlist.push_back(T(to1d(i,j,N), to1d(i-1,j,N), -2.0/h*y));
                        tripletlist.push_back(T(to1d(i,j,N), to1d(i-2,j,N), 1/2.0/h*y));
                    }
                    else
                    {
                        tripletlist.push_back(T(to1d(i,j,N), to1d(i+1,j,N), 1/2.0/h*y));
                        tripletlist.push_back(T(to1d(i,j,N), to1d(i-1,j,N), -1/2.0/h*y));
                    }

                    if(j == 0) // x direction
                    {
                        center_coeff += -3/2.0/h*x;
                        tripletlist.push_back(T(to1d(i,j,N), to1d(i,j+1,N), 2.0/h*x));
                        tripletlist.push_back(T(to1d(i,j,N), to1d(i,j+2,N), -1/2.0/h*x));
                    }
                    else if(j == N)
                    {
                        center_coeff += 3/2.0/h*x;
                        tripletlist.push_back(T(to1d(i,j,N), to1d(i,j-1,N), -2.0/h*x));
                        tripletlist.push_back(T(to1d(i,j,N), to1d(i,j-2,N), 1/2.0/h*x));
                    }
                    else
                    {
                        tripletlist.push_back(T(to1d(i,j,N), to1d(i,j+1,N), 1/2.0/h*x));
                        tripletlist.push_back(T(to1d(i,j,N), to1d(i,j-1,N), -1/2.0/h*x));
                    }
                    tripletlist.push_back(T(to1d(i,j,N), to1d(i,j,N), center_coeff));
                }
                else
                {
                    tripletlist.push_back(T(to1d(i,j,N), to1d(i,j,N), 4/SQUARE(h)));
                    tripletlist.push_back(T(to1d(i,j,N), to1d(i+1,j,N), -1/SQUARE(h)));
                    tripletlist.push_back(T(to1d(i,j,N), to1d(i-1,j,N), -1/SQUARE(h)));
                    tripletlist.push_back(T(to1d(i,j,N), to1d(i,j+1,N), -1/SQUARE(h)));
                    tripletlist.push_back(T(to1d(i,j,N), to1d(i,j-1,N), -1/SQUARE(h)));
                    if(N == originN)
                        _f(to1d(i,j,N)) = cond.f_val[i][j];
                }
            }
        A.setFromTriplets(tripletlist.begin(), tripletlist.end());
        _Matrix.push_back(A);
        N /= 2;
    }
    return 0;
}

int Generator<2>::gen_irregular(std::vector<SparseMatrix>& _Matrix, VectorXd &_f)
{
    int originN = cond.N;
    int N = originN;
    int s = 1;
    _f.resize(SQUARE(N+1));

    while(N >= coarsest_N)
    {
        SparseMatrix A(SQUARE(N+1), SQUARE(N+1));
        vecT tripletlist;
        double h = 1/double(N);

        for(int i = 0 ; i <= N ; i ++)
            for(int j = 0 ; j <= N ; j ++)
            {
                // std::cout << i << " " << j << std::endl;
                double rhs = 0;
                if(cond.loc[i*s][j*s] == 0) // out domain
                {
                    tripletlist.push_back(T(to1d(i,j,N), to1d(i,j,N), 1));
                    rhs = 0;
                }
                else if(cond.loc[i*s][j*s] == 1) // on boundary
                {
                    auto bv = cond.find_BV(Vec(j*h, i*h));
                    tripletlist.push_back(T(to1d(i,j,N), to1d(i,j,N), 1));
                    rhs = bv.val;
                }
                else // in domain
                {
                    double center_coeff = 0, tmp;
                    rhs = cond.f_val[i*s][j*s];
                    // left and right
                    double theta_left = 1, theta_right = 1, left, right;
                    if(cond.loc[i*s][(j+1)*s] == 0) // right is out domain
                    {
                        auto bv = cond.find_BV(Vec(j*h,i*h), Vec((j+1)*h,i*h));
                        theta_right = (bv.p[0] - j*h)/h;
                        right = bv.val;
                    }
                    if(cond.loc[i*s][(j-1)*s] == 0) // left is out domain
                    {
                        auto bv = cond.find_BV(Vec(j*h,i*h), Vec((j-1)*h,i*h));
                        theta_left = (j*h - bv.p[0])/h;
                        left = bv.val;
                    }
                    tmp = theta_left*theta_right*(theta_left+theta_right)/2.0*SQUARE(h);
                    center_coeff += (theta_left+theta_right)/tmp;
                    if(theta_left == 1) // left in domain
                        tripletlist.push_back(T(to1d(i,j,N), to1d(i,j-1,N),-theta_right/tmp));
                    else
                        rhs += theta_right/tmp*left;
                    if(theta_right == 1) // right in domain
                        tripletlist.push_back(T(to1d(i,j,N), to1d(i,j+1,N),-theta_left/tmp));
                    else
                        rhs += theta_left/tmp*right;
                    
                    // up and down
                    double theta_up = 1, theta_down = 1, up, down;
                    if(cond.loc[(i+1)*s][j*s] == 0) // up is out domain
                    {
                        auto bv = cond.find_BV(Vec(j*h,i*h), Vec(j*h,(i+1)*h));
                        theta_up = (bv.p[1] - i*h)/h;
                        up = bv.val;
                    }
                    // std::cout << "**" << std::endl;
                    if(cond.loc[(i-1)*s][j*s] == 0) // down is out domain
                    {
                        auto bv = cond.find_BV(Vec(j*h,i*h), Vec(j*h,(i-1)*h));
                        // std::cout << bv.p[0] << " " << bv.p[1] << std::endl;
                        theta_down = (i*h - bv.p[1])/h;
                        down = bv.val;
                    }
                    tmp = theta_up*theta_down*(theta_up+theta_down)/2.0*SQUARE(h);
                    center_coeff += (theta_up+theta_down)/tmp;
                    if(theta_down == 1) // down in domain
                        tripletlist.push_back(T(to1d(i,j,N), to1d(i-1,j,N),-theta_up/tmp));
                    else
                        rhs += theta_up/tmp*down;
                    if(theta_up == 1) // up in domain
                        tripletlist.push_back(T(to1d(i,j,N), to1d(i+1,j,N),-theta_down/tmp));
                    else
                        rhs += theta_down/tmp*up;
                    tripletlist.push_back(T(to1d(i,j,N), to1d(i,j,N),center_coeff));
                }
                
                if(N == originN)
                    _f(to1d(i,j,N)) = rhs;
            }
        A.setFromTriplets(tripletlist.begin(), tripletlist.end());
        _Matrix.push_back(A);
        N /= 2;
        s *= 2;
    }
    return 0;
}

/*
template<int dim>
int Generator<dim>::generate(std::vector<SparseMatrix>& _Matrix, VectorXd &_f)
{
    std::vector<SparseMatrix>().swap(_Matrix); //< clear the _Matrix
    if(dim == 1)
    {
        gen_dim1(_Matrix, _f);
    }
    else if(dim == 2)
    {
        if(cond.isRegular)
            gen_dim2_regular(_Matrix, _f);
        else
            gen_dim2_irregular(_Matrix, _f);
    }
    return 0;
}

template<int dim>
int Generator<dim>::gen_dim1(std::vector<SparseMatrix>& _Matrix, VectorXd &_f)
{
    int N = cond.N;
    _f.resize(N+1);
    SparseMatrix A(N+1, N+1);
    vecT tripletlist;
    for(int i = 0 ; i <= N ; i ++) // generate vector f
    {
        if(i == 0)
            _f(i) = cond.BV[0][0];
        else if(i == N)
            _f(i) = cond.BV[1][0];
        else
            _f(i) = cond.f_val[i];
    }
    while(N >= coarsest_N)
    {
        double h = 1/double(N);
        for(int i = 0 ; i <= N ; i ++)
        {
            if(i == 0)
            {
                if(cond.BV[0][1] == 0)
                {
                    tripletlist.push_back(T(0, 0, 1));
                }
                else // Neumann
                {
                    tripletlist.push_back(T(0, 0, -3/2.0/h));
                    tripletlist.push_back(T(0, 1, 2/h));
                    tripletlist.push_back(T(0, 2, -1/2.0/h));
                }
            }
            else if(i == N)
            {
                if(cond.BV[1][1] == 0)
                    tripletlist.push_back(T(N, N, 1));
                else // Neumann
                {
                    tripletlist.push_back(T(N, N, 3/2.0/h));
                    tripletlist.push_back(T(N, N-1, -2/h));
                    tripletlist.push_back(T(N, N-2, 1/2.0/h));
                }
            }
            else
            {
                tripletlist.push_back(T(i,i-1,1/SQUARE(h)));
                tripletlist.push_back(T(i,i,-2/SQUARE(h)));
                tripletlist.push_back(T(i,i+1,1/SQUARE(h)));
            }
        }
        A.setFromTriplets(tripletlist.begin(), tripletlist.end());
        _Matrix.push_back(A);
        N /= 2;
    }
    return 0;
}

template<int dim>
int Generator<dim>::gen_dim2_regular(std::vector<SparseMatrix>& _Matrix, VectorXd &_f)
{
    int originN = cond.N;
    int N = originN;
    _f.resize(SQUARE(N+1));
    
    while(N >= coarsest_N)
    {
        SparseMatrix A(SQUARE(N+1), SQUARE(N+1));
        vecT tripletlist;
        double h = 1/double(N);

        for(int i = 0 ; i <= N ; i ++)
            for(int j = 0 ; j <= N ; j ++)
            {
                if(i == 0 || j == 0 || i == N || j == N)
                {
                    auto bv = cond.find_BV(Vec(j*h, i*h));
                    if(N == originN)
                        _f(to1d(i,j,N)) = bv.val;
                    if(bv.n_vec == Vec(0,0)) // Dirichlet
                    {
                        tripletlist.push_back(T(to1d(i,j,N), to1d(i,j,N), 1));
                        continue;
                    }
                    // Neumann
                    double center_coeff = 0;
                    double x = bv.n_vec[0], y = bv.n_vec[1];
                    if(i == 0) // y direction
                    {
                        center_coeff += -3/2.0/h*y;
                        tripletlist.push_back(T(to1d(i,j,N), to1d(i+1,j,N), 2.0/h*y));
                        tripletlist.push_back(T(to1d(i,j,N), to1d(i+2,j,N), -1/2.0/h*y));
                    }
                    else if(i == N)
                    {
                        center_coeff += 3/2.0/h*y;
                        tripletlist.push_back(T(to1d(i,j,N), to1d(i-1,j,N), -2.0/h*y));
                        tripletlist.push_back(T(to1d(i,j,N), to1d(i-2,j,N), 1/2.0/h*y));
                    }
                    else
                    {
                        tripletlist.push_back(T(to1d(i,j,N), to1d(i+1,j,N), 1/2.0/h*y));
                        tripletlist.push_back(T(to1d(i,j,N), to1d(i-1,j,N), -1/2.0/h*y));
                    }

                    if(j == 0) // x direction
                    {
                        center_coeff += -3/2.0/h*x;
                        tripletlist.push_back(T(to1d(i,j,N), to1d(i,j+1,N), 2.0/h*x));
                        tripletlist.push_back(T(to1d(i,j,N), to1d(i,j+2,N), -1/2.0/h*x));
                    }
                    else if(j == N)
                    {
                        center_coeff += 3/2.0/h*x;
                        tripletlist.push_back(T(to1d(i,j,N), to1d(i,j-1,N), -2.0/h*x));
                        tripletlist.push_back(T(to1d(i,j,N), to1d(i,j-2,N), 1/2.0/h*x));
                    }
                    else
                    {
                        tripletlist.push_back(T(to1d(i,j,N), to1d(i,j+1,N), 1/2.0/h*x));
                        tripletlist.push_back(T(to1d(i,j,N), to1d(i,j-1,N), -1/2.0/h*x));
                    }
                    tripletlist.push_back(T(to1d(i,j,N), to1d(i,j,N), center_coeff));
                }
                else
                {
                    tripletlist.push_back(T(to1d(i,j,N), to1d(i,j,N), 4/SQUARE(h)));
                    tripletlist.push_back(T(to1d(i,j,N), to1d(i+1,j,N), -1/SQUARE(h)));
                    tripletlist.push_back(T(to1d(i,j,N), to1d(i-1,j,N), -1/SQUARE(h)));
                    tripletlist.push_back(T(to1d(i,j,N), to1d(i,j+1,N), -1/SQUARE(h)));
                    tripletlist.push_back(T(to1d(i,j,N), to1d(i,j-1,N), -1/SQUARE(h)));
                    if(N == originN)
                        _f(to1d(i,j,N)) = cond.f_val[i][j];
                }
            }
        A.setFromTriplets(tripletlist.begin(), tripletlist.end());
        _Matrix.push_back(A);
        N /= 2;
    }
    return 0;
}

template<int dim>
int Generator<dim>::gen_dim2_irregular(std::vector<SparseMatrix>& _Matrix, VectorXd &_f)
{
    int originN = cond.N;
    int N = originN;
    int s = 1;
    _f.resize(SQUARE(N+1));

    while(N >= coarsest_N)
    {
        SparseMatrix A(SQUARE(N+1), SQUARE(N+1));
        vecT tripletlist;
        double h = 1/double(N);

        for(int i = 0 ; i <= N ; i ++)
            for(int j = 0 ; j <= N ; j ++)
            {
                double rhs = 0;
                if(cond.loc[i*s][j*s] == 0) // out domain
                {
                    tripletlist.push_back(T(to1d(i,j,N), to1d(i,j,N), 1));
                    rhs = 0;
                }
                else if(cond.loc[i*s][j*s] == 1) // on boundary
                {
                    auto bv = cond.find_BV(Vec(j*h, i*h));
                    tripletlist.push_back(T(to1d(i,j,N), to1d(i,j,N), 1));
                    rhs = bv.val;
                }
                else // in domain
                {
                    double center_coeff = 0, tmp;
                    rhs = cond.f_val[i*s][j*s];
                    // left and right
                    double theta_left = 1, theta_right = 1, left, right;
                    if(cond.loc[i*s][(j+1)*s] == 0) // right is out domain
                    {
                        auto bv = cond.find_BV(Vec(j*h,i*h), Vec((j+1)*h,i*h));
                        theta_right = (bv.p[0] - j*h)/h;
                        right = bv.val;
                    }
                    if(cond.loc[i*s][(j-1)*s] == 0) // left is out domain
                    {
                        auto bv = cond.find_BV(Vec(j*h,i*h), Vec((j-1)*h,i*h));
                        theta_left = (j*h - bv.p[0])/h;
                        right = bv.val;
                    }
                    tmp = theta_left*theta_right*(theta_left+theta_right)/2.0*SQUARE(h);
                    center_coeff += (theta_left+theta_right)/tmp;
                    if(theta_left == 1) // left in domain
                        tripletlist.push_back(T(to1d(i,j,N), to1d(i,j-1,N),-theta_right/tmp));
                    else
                        rhs += theta_right/tmp*left;
                    if(theta_right == 1) // right in domain
                        tripletlist.push_back(T(to1d(i,j,N), to1d(i,j+1,N),-theta_left/tmp));
                    else
                        rhs += theta_left/tmp*right;

                    // up and down
                    double theta_up = 1, theta_down = 1, up, down;
                    if(cond.loc[(i+1)*s][j*s] == 0) // up is out domain
                    {
                        auto bv = cond.find_BV(Vec(j*h,i*h), Vec(j*h,(i+1)*h));
                        theta_up = (bv.p[1] - i*h)/h;
                        up = bv.val;
                    }
                    if(cond.loc[(i-1)*s][j*s] == 0) // down is out domain
                    {
                        auto bv = cond.find_BV(Vec(j*h,i*h), Vec(j*h,(i-1)*h));
                        theta_left = (i*h - bv.p[1])/h;
                        down = bv.val;
                    }
                    tmp = theta_up*theta_down*(theta_up+theta_down)/2.0*SQUARE(h);
                    center_coeff += (theta_up+theta_down)/tmp;
                    if(theta_down == 1) // down in domain
                        tripletlist.push_back(T(to1d(i,j,N), to1d(i-1,j,N),-theta_up/tmp));
                    else
                        rhs += theta_up/tmp*down;
                    if(theta_up == 1) // up in domain
                        tripletlist.push_back(T(to1d(i,j,N), to1d(i+1,j,N),-theta_down/tmp));
                    else
                        rhs += theta_down/tmp*up;
                }
                
                if(N == originN)
                    _f(to1d(i,j,N)) = rhs;
            }
        A.setFromTriplets(tripletlist.begin(), tripletlist.end());
        _Matrix.push_back(A);
        N /= 2;
        s *= 2;
    }
}
*/
// template class Generator<1>;
// template class Generator<2>;