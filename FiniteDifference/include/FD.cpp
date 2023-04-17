/**
 * @file FD.cpp
 * @author MaCheng (Ch.Ma01@outlook.com)
 * @brief implement of FD.h
 * @version 0.1
 * @date 2022-03-03
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include "FD.h"
#include <Eigen/Sparse>
#include <iostream>
#include <fstream>
#include "global.h"

// #define to1d(i,j) ((i)*(N+1)+(j))
// #define to2d1(n) (int(n/N))
// #define to2d2(n) (n-int(n/N)*N)

// #define isInCircle(p) (dist(p,center)<radius-EPS)

FDsolver::FDsolver(Conditions _cond): cond(_cond)
{
    int N = _cond.get_N();
    res.resize((N+1)*(N+1));
    if(!check_cond())
    {
        std::cerr << "invalid condition." << std::endl;
        return;
    }
}

bool FDsolver::check_cond()
{
    return true;
}

int FDsolver::general_solver()
{
    typedef Eigen::Triplet<double> T;

    const int N = cond.get_N();
    const double h = cond.get_h();
    auto center = cond.center;
    auto radius = cond.radius;
    //std::cout << N << " " << h << std::endl;
    Eigen::VectorXd b((N+1)*(N+1));
    Eigen::SparseMatrix<double> A((N+1)*(N+1), (N+1)*(N+1));
    std::vector<T> tripletlist;

    // for(auto it = (cond.BV0).begin(); it != (cond.BV0).end() ; it++)
    // {
    //     std::cout << it->first << " " << it->second << std::endl;
    // }
    // std::cout <<"***" <<  (cond.BV0).count(Vec(0.5, 1)) << std::endl;
    for(int i = 0 ; i <= N ; i ++)
        for(int j = 0 ; j <= N ; j ++)
        {
            Vec key(j*h, i*h);
            
            if(dist(key, center) < radius - EPS) // in circle
            {
                cond.f_val[i][j] = D_MAX;
                tripletlist.push_back(T(to1d(i,j),to1d(i,j),1));
                b(to1d(i,j)) = 0;
                continue;
            }
            // if(cond.f_val[i][j] > 10000)
            // {
            //     std::cout << "err: " << key << std::endl;
            // }
            double delta_RHS = 0;
            std::vector<double> coeff((N+1)*(N+1), 0);

            // boundary of square
            if(i == 0 || j == 0 || i == N || j == N)
            {
                auto bv = cond.find_BV(key);
                if(bv.n_vec == Vec(0,0)) // value is given
                {
                    // std::cout << i << " " << j << std::endl;
                    tripletlist.push_back(T(to1d(i,j),to1d(i,j),1));
                    b(to1d(i,j)) = bv.val;
                }
                else
                {
                    // std::cout << i << " " << j << std::endl;
                    if(isInSquare(RIGHT(key)))
                    {
                        if(isInSquare(LEFT(key)))
                        {
                            if(dist(LEFT(key),center) > radius-EPS) // left out circle
                            {
                                if(dist(RIGHT(key), center) > radius-EPS) // both out circle
                                {
                                    // std::cout << "LR" << std::endl;
                                    coeff[to1d(i,j+1)] += bv.n_vec[0]/2/h;
                                    coeff[to1d(i,j-1)] += -bv.n_vec[0]/2/h;
                                }
                                else // only left out circle
                                {
                                    coeff[to1d(i,j-1)] += -bv.n_vec[0]/h;
                                    coeff[to1d(i,j)] += bv.n_vec[0]/h;
                                }
                            }
                            else // only right out circle
                            {
                                coeff[to1d(i,j)] += -bv.n_vec[0]/h;
                                coeff[to1d(i,j+1)] += bv.n_vec[0]/h;
                            }
                        }
                        else // only right in square
                        {
                            if(dist(RIGHT(key), center) > radius-EPS) // right out circle
                            {
                                // std::cout << "R" << std::endl;
                                coeff[to1d(i,j)] += -3*bv.n_vec[0]/2/h;
                                coeff[to1d(i,j+1)] += 2*bv.n_vec[0]/h;
                                coeff[to1d(i,j+2)] += -1*bv.n_vec[0]/2/h;
                            }
                            else // right irregular
                            {
                                double t = center[0] - std::sqrt(SQUARE(radius)-SQUARE(i*h-center[1]));
                                auto bv1 = cond.find_BV(Vec(t, i*h));
                                double theta = ((bv1.p)[0] - key[0]) / h;
                                if(bv1.n_vec == Vec(0,0))
                                {
                                    coeff[to1d(i,j)] += -bv.n_vec[0]/(theta*h);
                                    delta_RHS += -bv.n_vec[0]*bv1.val/(theta*h);
                                }
                                else
                                {
                                    int i_star = cond.find_i(i, bv1.p[0]);
                                    int delta_i = i_star - i;
                                    double A1 = bv.n_vec[0]/(theta*h);
                                    double A2 = bv1.n_vec[1]/(delta_i*h) - bv1.n_vec[0]/(theta*h);
                                    coeff[to1d(i_star, j)] += (1-theta)/(delta_i*h)/A2*A1*bv1.n_vec[1];
                                    coeff[to1d(i_star, j+1)] += theta/(delta_i*h)/A2*A1*bv1.n_vec[1];
                                    coeff[to1d(i,j)] += -bv1.n_vec[0]/(theta*h)/A2*A1-A1;
                                    delta_RHS += bv1.val/A2*A1;
                                }   
                            }
                        }
                    }
                    else // right out square, then left must in square
                    {
                        if(dist(LEFT(key),center) > radius-EPS)
                        {
                            // std::cout << "L" << std::endl;
                            coeff[to1d(i,j-2)] += 1*bv.n_vec[0]/2/h;
                            coeff[to1d(i,j-1)] += -2*bv.n_vec[0]/h;
                            coeff[to1d(i,j)] += 3*bv.n_vec[0]/2/h;
                        }
                        else // left irregular
                        {
                            double t = center[0] + std::sqrt(SQUARE(radius)-SQUARE(i*h-center[1]));
                            auto bv1 = cond.find_BV(Vec(t, i*h));
                            double theta = (key[0] - (bv1.p)[0]) / h;
                            if(bv1.n_vec == Vec(0,0))
                            {
                                coeff[to1d(i,j)] += bv.n_vec[0]/(theta*h);
                                delta_RHS += bv.n_vec[0]*bv1.val/(theta*h);
                            }
                            else
                            {
                                int i_star = cond.find_i(i, bv1.p[0]);
                                int delta_i = i_star - i;
                                double A1 = bv.n_vec[0]/(theta*h);
                                double A2 = bv1.n_vec[1]/(delta_i*h) + bv1.n_vec[0]/(theta*h);
                                coeff[to1d(i_star, j)] += -(1-theta)/(delta_i*h)/A2*A1*bv1.n_vec[1];
                                coeff[to1d(i_star, j+1)] += -theta/(delta_i*h)/A2*A1*bv1.n_vec[1];
                                coeff[to1d(i,j)] += -bv1.n_vec[0]/(theta*h)/A2*A1+A1;
                                delta_RHS += -bv1.val/A2*A1;
                            }
                        }
                        
                    }
                    
                    if(isInSquare(UP(key))) // UP in square
                    {
                        if(isInSquare(DOWN(key)))
                        {
                            if(dist(DOWN(key),center) > radius-EPS) // down out circle
                            {
                                if(dist(UP(key), center) > radius-EPS) // both out circle
                                {
                                    // std::cout << "UD" << std::endl;
                                    coeff[to1d(i+1,j)] += bv.n_vec[1]/2/h;
                                    coeff[to1d(i-1,j)] += -bv.n_vec[1]/2/h;
                                }
                                else // only down out circle
                                {
                                    coeff[to1d(i,j)] += bv.n_vec[1]/h;
                                    coeff[to1d(i-1,j)] += -bv.n_vec[1]/h;
                                }
                            }
                            else // only up out circle
                            {
                                coeff[to1d(i+1,j)] += bv.n_vec[1]/h;
                                coeff[to1d(i,j)] += -bv.n_vec[1]/h;
                            }
                        }
                        else // only up in square
                        {
                            if(dist(UP(key),center)>radius-EPS)
                            {
                                // std::cout << "U" << std::endl;
                                coeff[to1d(i+2,j)] += -1*bv.n_vec[1]/2/h;
                                coeff[to1d(i+1,j)] += 2*bv.n_vec[1]/h;
                                coeff[to1d(i,j)] += -3*bv.n_vec[1]/2/h;
                            }
                            else // irregular up
                            {
                                double t = center[1] - std::sqrt(SQUARE(radius)-SQUARE(j*h-center[0]));
                                auto bv1 = cond.find_BV(Vec(j*h, t));
                                double theta = ((bv1.p)[1] - key[1]) / h;
                                if(bv1.n_vec == Vec(0,0)) // value is given
                                {
                                    coeff[to1d(i,j)] += -bv.n_vec[1]/(theta*h);
                                    delta_RHS += -bv1.val*bv.n_vec[1]/(theta*h);
                                } 
                                else
                                {
                                    int j_star = cond.find_j(j, bv.p[1]);
                                    int delta_j = j_star - j; 
                                    double A1 = bv.n_vec[1]/(theta*h);
                                    double A2 = -bv1.n_vec[1]/(theta*h) + bv1.n_vec[0]/(delta_j*h);
                                    coeff[to1d(i, j_star)] += (1-theta)*bv1.n_vec[0]/(delta_j*h)/A2*A1;
                                    coeff[to1d(i+1, j_star)] += theta*bv1.n_vec[0]/(delta_j*h)/A2*A1;
                                    coeff[to1d(i,j)] += -bv1.n_vec[1]/(theta*h)/A2*A1-A1;
                                    delta_RHS += bv1.val/A2*A1;
                                }
                            }
                        }
                    }
                    else // only down in square
                    {
                        if(dist(DOWN(key),center)>radius-EPS)
                        {
                            // std::cout << "D" << std::endl;
                            coeff[to1d(i,j)] += 3*bv.n_vec[1]/2/h;
                            coeff[to1d(i-1,j)] += -2*bv.n_vec[1]/h;
                            coeff[to1d(i-2,j)] += 1*bv.n_vec[1]/2/h;
                        }
                        else // irregular down
                        {
                            double t = center[1] + std::sqrt(SQUARE(radius)-SQUARE(j*h-center[0]));
                            auto bv1 = cond.find_BV(Vec(j*h, t));
                            double theta = (key[1] - (bv1.p)[1]) / h;
                            if(bv1.n_vec == Vec(0,0)) // value is given
                            {
                                coeff[to1d(i,j)] += bv.n_vec[1]/(theta*h);
                                delta_RHS += bv1.val*bv.n_vec[1]/(theta*h);
                            } 
                            else
                            {
                                int j_star = cond.find_j(j, bv.p[1]);
                                int delta_j = j_star - j; 
                                double A1 = bv.n_vec[1]/(theta*h);
                                double A2 = bv1.n_vec[1]/(theta*h) + bv1.n_vec[0]/(delta_j*h);
                                coeff[to1d(i, j_star)] += -(1-theta)*bv1.n_vec[0]/(delta_j*h)/A2*A1;
                                coeff[to1d(i+1, j_star)] += -theta*bv1.n_vec[0]/(delta_j*h)/A2*A1;
                                coeff[to1d(i,j)] += -bv1.n_vec[1]/(theta*h)/A2*A1+A1;
                                delta_RHS += -bv.val/A2*A1;
                            }
                        }
                        
                    }
                    for(int k = 0 ; k < SQUARE(N+1) ; k ++)
                        if(coeff[k] != 0)
                            tripletlist.push_back(T(to1d(i,j),k,coeff[k]));
                    b(to1d(i,j)) = bv.val + delta_RHS; 
                }
                continue;
            }
            
            // on the circle
            if(std::abs(dist(key, center) - radius) < EPS)
            {
                auto bv = cond.find_BV(key);
                if(bv.n_vec == Vec(0,0)) // value is given
                {
                    tripletlist.push_back(T(to1d(i,j),to1d(i,j),1));
                    b(to1d(i,j)) = bv.val;
                }
                else
                {
                    double coeff_center = 0;
                    if(dist(center, RIGHT(key)) > radius-EPS) // right point is outside 
                    {
                        tripletlist.push_back(T(to1d(i,j),to1d(i,j+1),bv.n_vec[0]/h)); // right
                        coeff_center += -bv.n_vec[0]/h;
                    }
                    else
                    {
                        tripletlist.push_back(T(to1d(i,j),to1d(i,j-1),-bv.n_vec[0]/h)); // left
                        coeff_center += bv.n_vec[0]/h;
                    }
                    if(dist(center, DOWN(key)) > radius-EPS)
                    {
                        tripletlist.push_back(T(to1d(i,j),to1d(i-1,j),-bv.n_vec[1]/h)); // down
                        coeff_center += bv.n_vec[1]/h;
                    }
                    else
                    {
                        tripletlist.push_back(T(to1d(i,j),to1d(i+1,j),bv.n_vec[1]/h)); // up
                        coeff_center += -bv.n_vec[0]/h;
                    }
                    tripletlist.push_back(T(to1d(i,j),to1d(i,j),coeff_center));
                    b(to1d(i,j)) = bv.val;
                }
                continue;
            }
            
            
            // inner point
            if(dist(center, LEFT(key)) < radius-EPS) // left point isn't in domain.
            {
                ////std::cout<<"ir"<<std::endl;
                double t = center[0] + std::sqrt(SQUARE(radius)-SQUARE(i*h-center[1]));
                auto bv = cond.find_BV(Vec(t, i*h));
                double theta = (key[0] - (bv.p)[0]) / h;
                double A1 = 0.5*theta*(1+theta)*SQUARE(h);
                if(bv.n_vec == Vec(0,0)) // value is given
                {
                    coeff[to1d(i,j)] += (1+theta)/A1;
                    coeff[to1d(i,j+1)] += -theta/A1;
                    delta_RHS += bv.val/A1; 
                }
                else
                {
                    ////std::cout<<"ir"<<std::endl;
                    int i_star = cond.find_i(i, bv.p[0]);
                    // std::cout << i_star << std::endl;
                    int delta_i = i_star - i;
                    double A2 = bv.n_vec[1]/(delta_i*h) + bv.n_vec[0]/(theta*h);
                    coeff[to1d(i_star, j)] += -(1-theta)/(delta_i*h)/A2/A1*bv.n_vec[1];
                    coeff[to1d(i_star, j-1)] += -theta/(delta_i*h)/A2/A1*bv.n_vec[1];
                    coeff[to1d(i,j)] += -bv.n_vec[0]/(theta*h)/A2/A1+(1+theta)/A1;
                    coeff[to1d(i,j+1)] += -theta/A1;
                    delta_RHS += -bv.val/A2/A1;
                }
            }
            else if(dist(center, RIGHT(key)) < radius-EPS) // right point isn't in domain
            {
                ////std::cout<<"ir"<<std::endl;
                double t = center[0] - std::sqrt(SQUARE(radius)-SQUARE(i*h-center[1]));
                auto bv = cond.find_BV(Vec(t, i*h));
                double theta = ((bv.p)[0] - key[0]) / h;
                double A1 = 0.5*theta*(1+theta)*SQUARE(h);
                if(bv.n_vec == Vec(0,0)) // value is given
                {
                    coeff[to1d(i,j)] += (1+theta)/A1;
                    coeff[to1d(i,j-1)] += -theta/A1;
                    delta_RHS += bv.val/A1; 
                } 
                else
                {
                    ////std::cout<<"ir"<<std::endl;
                    int i_star = cond.find_i(i, bv.p[0]);
                    int delta_i = i_star - i;
                    double A2 = bv.n_vec[1]/(delta_i*h) - bv.n_vec[0]/(theta*h);
                    coeff[to1d(i_star, j)] += -(1-theta)/(delta_i*h)/A2/A1*bv.n_vec[1];
                    coeff[to1d(i_star, j+1)] += -theta/(delta_i*h)/A2/A1*bv.n_vec[1];
                    coeff[to1d(i,j)] += bv.n_vec[0]/(theta*h)/A2/A1+(1+theta)/A1;
                    coeff[to1d(i,j-1)] += -theta/A1;
                    delta_RHS += -bv.val/A2/A1;
                }
            }
            else // both left and right points are in domain.
            {
                coeff[to1d(i,j+1)] += -1/SQUARE(h);
                coeff[to1d(i,j-1)] += -1/SQUARE(h);
                coeff[to1d(i,j)] += 2/SQUARE(h);
            }
            
            if(dist(center, UP(key)) < radius-EPS) // up point isn't in domain
            {
                ////std::cout<<"ir"<<std::endl;
                double t = center[1] - std::sqrt(SQUARE(radius)-SQUARE(j*h-center[0]));
                auto bv = cond.find_BV(Vec(j*h, t));
                double theta = ((bv.p)[1] - key[1]) / h;
                double A1 = 0.5*theta*(1+theta)*SQUARE(h);
                if(bv.n_vec == Vec(0,0)) // value is given
                {
                    coeff[to1d(i,j)] += (1+theta)/A1;
                    coeff[to1d(i-1,j)] += -theta/A1;
                    delta_RHS += bv.val/A1; 
                } 
                else
                {
                    //std::cout<<"ir"<<std::endl;
                    int j_star = cond.find_j(j, bv.p[1]);
                    int delta_j = j_star - j; 
                    double A2 = -bv.n_vec[1]/(theta*h) + bv.n_vec[0]/(delta_j*h);
                    coeff[to1d(i, j_star)] += -(1-theta)*bv.n_vec[0]/(delta_j*h)/A2/A1;
                    coeff[to1d(i+1, j_star)] += -theta*bv.n_vec[0]/(delta_j*h)/A2/A1;
                    coeff[to1d(i,j)] += bv.n_vec[1]/(theta*h)/A2/A1+(1+theta)/A1;
                    coeff[to1d(i-1,j)] += -theta/A1;
                    delta_RHS += -bv.val/A2/A1;
                }
            }
            else if(dist(center, DOWN(key)) < radius-EPS) // down point isn't in domain
            {
                ////std::cout<<"ir"<<std::endl;
                double t = center[1] + std::sqrt(SQUARE(radius)-SQUARE(j*h-center[0]));
                auto bv = cond.find_BV(Vec(j*h, t));
                double theta = (key[1] - (bv.p)[1]) / h;
                double A1 = 0.5*theta*(1+theta)*SQUARE(h);
                if(bv.n_vec == Vec(0,0)) // value is given
                {
                    coeff[to1d(i,j)] += (1+theta)/A1;
                    coeff[to1d(i+1,j)] += -theta/A1;
                    delta_RHS += bv.val/A1; 
                } 
                else
                {
                    int j_star = cond.find_j(j, bv.p[1]);
                    int delta_j = j_star - j;
                    double A2 = bv.n_vec[1]/(theta*h) + bv.n_vec[0]/(delta_j*h);
                    coeff[to1d(i, j_star)] += -(1-theta)*bv.n_vec[0]/(delta_j*h)/A2/A1;
                    coeff[to1d(i-1, j_star)] += -theta*bv.n_vec[0]/(delta_j*h)/A2/A1;
                    coeff[to1d(i,j)] += -bv.n_vec[1]/(theta*h)/A2/A1+(1+theta)/A1;
                    coeff[to1d(i+1,j)] += -theta/A1;
                    delta_RHS += -bv.val/A2/A1;
                    // if(j==2 && i==3)
                    // {
                    //     std::cout << A1 << " " << A2 << " " << delta_RHS << std::endl;
                    // }
                }
            }
            else // both up and down points are in domain.
            {
                coeff[to1d(i+1,j)] += -1/SQUARE(h);
                coeff[to1d(i-1,j)] += -1/SQUARE(h);
                coeff[to1d(i,j)] += 2/SQUARE(h);
            }/**/
            // std::cout << i << " " << j << " " << cond.f_val[i][j] << " " << delta_RHS << std::endl;
            for(int k = 0 ; k < SQUARE(N+1) ; k ++)
                if(coeff[k] != 0)
                    tripletlist.push_back(T(to1d(i,j),k,coeff[k]));
            b(to1d(i,j)) = cond.f_val[i][j] + delta_RHS; 
        }
    A.setFromTriplets(tripletlist.begin(), tripletlist.end());
    // solve the linear system
    // for(int i = 0 ; i < SQUARE(N+1) ; i ++)
    // {
    //     for(int j = 0 ; j < SQUARE(N+1) ; j ++)
    //         std::cout << A(i,j) << " ";
    //     std::cout << std::endl;
    // // }
    // if(output)
    // {
    //     std::cout << A << std::endl;
    //     for(int i = 0 ; i < SQUARE(N+1) ; i ++)
    //         std::cout << b(i) << std::endl;
    // }
    
    // for(int i = 0 ; i <= N+1 ; i ++)
    //     for(int j = 0 ; j <= N + 1 ; j ++)
    //         std::cout <<i << " " << j << " " <<  b(to1d(i,j)) << std::endl;
    
    A.makeCompressed();
    Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > solver;
    solver.analyzePattern(A);
    solver.factorize(A);
    res = solver.solve(b);
    return 0;
}

int FDsolver::get_err(double(*u)(double, double), std::string _path, int _index = 0)
{
    int N = cond.get_N();
    double h = cond.get_h();
    double1d abs_err, rel_err;
    for(int i = 0; i <= N ; i ++)
        for(int j = 0 ; j <= N ; j ++)
        {
            if(dist(Vec(j*h, i*h), cond.center) < cond.radius - EPS)
            {
                abs_err.push_back(0);
                rel_err.push_back(0);
                continue;
            }
            double e = std::abs(res(to1d(i,j)) - u(j*h, i*h));
            abs_err.push_back(e);
            rel_err.push_back(e/std::abs(u(j*h, i*h)));
        }
    std::ofstream outfile;
    outfile.open(_path);
    if(!outfile.is_open()) 
    {
        std::cerr << "invalid path." << std::endl;
        return -1;
    }
    outfile << "x = [ ";
    for(int i = 0 ; i <= N ; i ++)
        for(int j = 0 ; j <= N ; j ++)
            outfile << j*h << " ";
    outfile << "];"<< std::endl;

    outfile << "y = [ ";
    for(int i = 0 ; i <= N ; i ++)
        for(int j = 0 ; j <= N ; j ++)
            outfile << i*h << " ";
    outfile << "];"<< std::endl;

    outfile << "abs_err = [ ";
    for(int i = 0; i <= N ; i ++)
        for(int j = 0 ; j <= N ; j ++)
            outfile << abs_err[to1d(i,j)] << " ";
    outfile << "];" << std::endl;

    outfile << "rel_err = [ ";
    for(int i = 0; i <= N ; i ++)
        for(int j = 0 ; j <= N ; j ++)
            outfile << rel_err[to1d(i,j)] << " ";
    outfile << "];" << std::endl;

    outfile << "scatter3(x, y, abs_err, 50, abs_err, '.');" << std::endl << "colorbar;" << std::endl;
    outfile << "title('Test "<< _index << ":$N = "<<N<<"$','Interpreter','latex','FontSize',20);" << std::endl;
    outfile << "saveas(gcf, './fig/f" << _index << "_N" << N << ".png');" << std::endl;
    //title('Test 1: $N = 64$','Interpreter','latex','FontSize',20);
    //saveas(gcf, './fig/f5N64.png');
    outfile.close();
    return 0;
}

double FDsolver::get_err(double(*u)(double, double), int _type)
{
    int N = cond.get_N();
    double h = cond.get_h();
    double e1 = 0 , e2 = 0 , einf = 0;
    for(int i = 0; i <= N ; i ++)
        for(int j = 0 ; j <= N ; j ++)
        {
            if(dist(Vec(j*h, i*h), cond.center) < cond.radius - EPS)
                continue;
            double err = std::abs(res(to1d(i,j)) - u(j*h, i*h));
            einf = max(err, einf);
            e1 += err;
            e2 += err*err;
        }
    e1 = e1 * h * h;
    e2 = std::sqrt(h*h*e2);
    switch(_type)
    {
        case 0: return einf;
            break;
        case 1: return e1;
            break;
        case 2: return e2;
            break;
    }
    return -1;
}

std::ostream &operator<<(std::ostream &os, const FDsolver &_p)
{
    int N = (_p.cond).get_N();
    int k = 0;
    for(int i = 0 ; i <= N ; i ++)
    {
        for(int j = 0 ; j <= N ; j ++)
            std::cout << _p.res[k++] << " ";
        std::cout << std::endl;
    }
    return os;
}