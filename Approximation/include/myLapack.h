/**
 * @file myLapack.h
 * @author Ouyang Shangke
 * @brief Funs in this file are mainly for solving linear system, by invoking the func in "lapacke.h".
 * @details The aim of those funcs is re-encapsulation the original function in "lapacke.h", 
 * so that one can use it conviently instead of declaring so many parameters before invoking it.
 * @version 0.1
 * @date 2021-11-23
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#ifndef __MYLAPCK_H__
#define __MYLAPCK_H__

#include "Config.h"
#include <lapacke.h>

/**
 * @brief Transform a vector to a simple array.
 * 
 * @tparam TYPE type of the vector
 * @param _t the vector
 * @return array
 */
template<class TYPE>
inline TYPE* toArray (std::vector<TYPE> _t)
{
    TYPE* res = new TYPE[_t.size()];
    int i = 0;
    for(auto it = _t.begin(); it < _t.end(); it ++, i ++)
        res[i] = *it;
    return res;
}

/**
 * @brief Solve a linear system whose coefficient matrix is tridiagonal.
 * the core of this func is by "LAPACKE_dgtsv" in "Lapcke.h"
 * @param _DL sub-diagonal of coefficient matrix
 * @param _D diagonal of coefficient matrix
 * @param _DU super-diagonal of coefficient matrix
 * @param _b(in & out) input as a column vector, 
 *          and it will store the solution to this linear system after invoking this func.
 * @return 0 for OK
 */
int solve_tridiag(vector_double &_DL,vector_double &_D,vector_double &_DU,vector_double &_b)
{
    lapack_int N = _D.size();
    lapack_int NRHS = 1;
    lapack_int LDB = _b.size();
    lapack_int INFO;
    double* DL = toArray(_DL);
    double* D = toArray(_D);
    double* DU = toArray(_DU);
    double* b = toArray(_b);

    INFO = LAPACKE_dgtsv(LAPACK_COL_MAJOR, N, NRHS, DL, D, DU, b, LDB);
    
    if(INFO != 0)
        return INFO;
    
    for(int i = 0; i < N ; i ++)
        _b[i] = b[i];
    
    delete[] DL, D, DU, b;
    return 0;
}

/**
 * @brief Solve a general linear system. 
 * the core of this func is by "LAPACKE_dgesv" in "Lapcke.h"
 * @param _A store the coefficient matrix with column-major order
 * @param _b (in & out) input as a column vector, 
 *           and it will store the solution to this linear system after invoking this func.
 * @return 0 for OK
 */
int solve_general(vector2_double &_A, vector_double &_b)
{
    lapack_int N = _A.size();
    lapack_int NRHS = 1;
    lapack_int LDA = _A.size();
    lapack_int LDb = 1;
    lapack_int* IPIV = new int[N];
    lapack_int INFO;
    double *b = toArray(_b);

    ///< Transform the coef matrix
    double *A = new double[N*N];
    auto it = _A.begin();
    for(int i = 0 ; i < N && it < _A.end() ; i ++, it++)
    {
        for(int j = 0 ; j < N ; j ++)
        {
            if(j < (*it).size())
                A[i*N + j] = (*it)[j];
            else
                A[i*N + j] = 0;
        }
    }
    INFO = LAPACKE_dgesv(LAPACK_ROW_MAJOR, N, NRHS, A, LDA, IPIV, b, LDb);

    if(INFO != 0)
        return INFO;
    
    for(int i = 0; i < N ; i ++)
        _b[i] = b[i];

    delete[] A, b, IPIV;
    return 0;
}

/**
 * @brief QR-factorization 
 * 
 * @param A Original matrix.
 * @param Q [out]
 * @param R [out]
 * @return 0 for OK
 */
int QR(vector2_double& A, vector2_double& Q, vector2_double& R)
{
    int n = A.size();
    int m = A[0].size();
    Q.resize(n); R.resize(n);
    for(int i = 0 ; i < n ; i ++)
    {
        Q[i].resize(n);
        R[i].resize(m);
    }
    for(int i = 0 ; i < m ; i ++)
    {
        double sum = 0;
        for(int j = 0 ; j < n ; j ++)
            sum += A[j][i] * A[j][i];
        R[i][i] = std::sqrt(sum);
        for(int j = 0 ; j < n ; j ++)
            Q[j][i] = A[j][i] / R[i][i];
        for(int j = i + 1 ; j < m ; j ++)
        {
            double t = 0;
            for(int k = 0 ; k < n ; k ++)
                t += A[k][j] * Q[k][i];
            R[i][j] = t;
            for(int k = 0 ; k < n ; k ++)
                A[k][j] = A[k][j] - R[i][j] * Q[k][i];
        }
    }
    return 0;
}


#endif