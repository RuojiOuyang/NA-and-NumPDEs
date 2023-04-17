/**
 * @file Polynomial.h
 * @author OuyangShangke
 * @brief This file is for the class Polynomial.
 * @version 0.1
 * @date 2021-11-23
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#ifndef __Polynomial_H__
#define __Polynomial_H__

#include "Config.h"

/**
 * @brief class Polynomial.
 * 
 * @tparam _Degree degree of polynomial
 * @tparam _Dim dimention of polynomial
 */
template<u_int _Degree, u_int _Dim>
class Polynomial
{
private:

    /**
     * CoefType: The type of the coefficients.
     * vector_Coef: A vector stores the coefficients.
     */
    typedef Vec<double, _Dim> CoefType;
    typedef std::vector<Vec<double,_Dim>> vector_Coef;

    std::vector<CoefType> Coef; ///< coefficients of polynomial

public:
    enum{
        Dim = _Dim, 
        Degree = _Degree 
    };

    /**
     * @brief Default constructor.
     * Construct a new Polynomial object, and resize the vector Coef to _Degree+1.
     */
    Polynomial();

    /**
     * @brief Construct a new Polynomial object by a set of given coefficients.
     * 
     * @param _Coef given coefficients.
     */
    Polynomial(const vector_Coef& _Coef);

    /**
     * @brief Construct a new Polynomial object by a given Polynomial object.
     * 
     * @param _p given Polynomial object
     */
    Polynomial(const Polynomial& _p);

    /**
     * @brief Set the Coef object by a given vector.
     * 
     * @param _Coef given vector
     * @return
     *  @retval 0 for OK
     *  @retval -1 for error(may due to the wrong size of _Coef)
     */
    int set_Coef(const vector_Coef& _Coef);

    /**
     * @brief Set the coefficient of the polynomial at a certain degree.
     * 
     * @param _Deg given Degree
     * @param _Coef given coef
     * @return
     *  @retval 0 for OK
     *  @retval -1 for error(may due to the _Degree out of range)
     */
    int set_Coef(int _Deg, const CoefType& _Coef);

    /**
     * @brief Get the coefficient at a certain degree.
     * 
     * @param _Deg given degree
     * @return the corresponding coefficient
     */
    CoefType get_Coef(int _Deg);

    /**
     * @brief Get the derivative of the polynomial.
     * 
     * @return the derivative
     */
    Polynomial<_Degree==0?0:_Degree-1, _Dim> diff();

    /**
     * @brief Add two Polynomial objects.
     * 
     * @tparam _Degree2 the Degree of the other polynomial.
     * @param _p the other polynomial
     * @return result(a polynomial whose Degree equals to the maximum of the 2 polynomials Degrees)
     */
    template<u_int _Degree2>
    Polynomial<std::max(_Degree, _Degree2), _Dim> operator+(const Polynomial<_Degree2, _Dim>& _p) const;

    /**
     * @brief Substraction of two polynomial.
     * 
     * @tparam _Degree2 the Degree of the other polynomial
     * @param _p the other polynomial
     * @return result(a polynomial whose Degree equals to the maximum of the 2 polynomials Degrees)
     */
    template<u_int _Degree2>
    Polynomial<std::max(_Degree, _Degree2), _Dim> operator-(const Polynomial<_Degree2, _Dim>& _p) const;

    /**
     * @brief Multiplication of two polynomial.
     * 
     * @tparam _Degree2 The Degree of the other polynomial.
     * @param _p the other polynomial
     * @return result(a polynomial whose Degree equals to the sum of the 2 polynomials Degrees)
     */
    template<u_int _Degree2>
    Polynomial<_Degree+_Degree2, _Dim> operator*(const Polynomial<_Degree2, _Dim>& _p) const;

    /**
     * @brief Get the value of the polynomial at given point.
     * 
     * @param _x given point
     * @return function value
     */
    CoefType operator()(double _x) const;

    /**
     * @brief Output a polynomial.
     * 
     * @param output 
     * @param _p the polynomial
     * @return
     */
    template<u_int __Degree, u_int __Dim>
    friend std::ostream& operator<< (std::ostream& output, const Polynomial<__Degree, __Dim>& _p);
};

/// IMPLEMENT
#include "../src/header/Polynomial.cpp"

#endif