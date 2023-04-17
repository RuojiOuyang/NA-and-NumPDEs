#ifndef __PolyInterp_H__
#define __PolyInterp_H__

#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <vector>
#include <limits>

#define double_array std::vector<double>
#define double2_array std::vector<std::vector<double>>
#define double_array_it std::vector<double>::iterator
#define double_array_const_it std::vector<double>::const_iterator
#define D_MAX (std::numeric_limits<double>::max())
#define int_array std::vector<int>

struct InterpConditions
{
    double_array x; ///< interpolation sites
    double2_array fx; ///< corresponding function values

    /**
     * @brief Get the maximum order of the derivatives.
     * 
     * @return maximum order of the derivatives
     */
    int get_max_diff() const;
    /**
     * @brief Get the amount of all the function and its derivatives values.
     * 
     * @return total amount
     */
    int get_N() const;
    /**
     * @brief Check if the given point is one of the interpolation points.
     * 
     * @param _x given point
     * @return 
     *  @retval >=0 its position
     *  @retval -1 its not a interpolation point
     */
    int find_x(double _x) const;
    /**
     * @brief Input a InterpCondtions by a given structure.
     * 
     * @param input 
     * @param _interp given structure
     * @return
     */
    friend std::istream& operator>> (std::istream& input, InterpConditions& _interp);
};

class Polynomial
{
private:
    double_array coefficient; ///< coefficients of polynomial

public:
    /**
     * @brief Default constructor.
     * Construct a new Polynomial object.
     */
    Polynomial();

    /**
     * @brief Construct a new Polynomial object by a set of given coefficients.
     * 
     * @param _coefficient given coefficients
     */
    Polynomial(const double_array& _coefficient);

    /**
     * @brief Construct a new Polynomial object by a given Polynomial object.
     * 
     * @param _p given Polynomial object
     */
    Polynomial(const Polynomial& _p);

    int set_coefficient(const double_array& _coefficient);

    /**
     * @brief Get the degree of the polynomial.
     * 
     * @return the degree
     */
    int get_degree() const;

    /**
     * @brief Delete front zeros in coefficient.
     * 
     * @return 0 for OK
     */
    int refresh();

    /**
     * @brief Get the derivative of the polynomial.
     * 
     * @return the derivative
     */
    Polynomial diff();

    /**
     * @brief Get the derivative of given order.
     * 
     * @param _n given order
     * @return the derivative of given order
     */
    Polynomial diff(int _n);

    /**
     * @brief Output the polynomial to a given file in certain format,
     * so that we can use the expression when coding in Matlab or Octave.
     * 
     * @param path the path of the given file
     * @param _n a subscript of the polynomial
     * @return 0 for OK
     */
    int output4Matlab(const std::string& path, int _n = 0);


    /**
     * @brief Assign the object by a given Polynomial object.
     * 
     * @param _p given Polynomial object.
     * @return itself
     */
    Polynomial& operator= (const Polynomial& _p);
    /**
     * @brief Assign the object by a given real number 
     * to get a constant polynomial.
     * 
     * @param c given real number
     * @return itselt
     */
    Polynomial& operator= (const double c);

    /**
     * @brief Determine whether two polynomial is equivalent.
     * 
     * @param _p the other polynomial
     * @return 
     *  @retval true for equal
     *  @retval false for unequal
     */
    bool operator== (const Polynomial& _p) const;

    /**
     * @brief Add two Polynomial objects.
     * 
     * @param _p the oter polynomial
     * @return result
     */
    Polynomial operator+ (const Polynomial& _p) const;
    /**
     * @brief Add the polynomial with a real number.
     * 
     * @param _t the given real number
     * @return result
     */
    Polynomial operator+ (double _t) const;
    /**
     * @brief Add a real number with a given polynomial
     * 
     * @param _t real number
     * @param _p given polynomial
     * @return result
     */
    friend Polynomial operator+ (double _t, const Polynomial& _p);

    /**
     * @brief Get the opposite polynomial of this polynomial.
     * 
     * @return result
     */
    Polynomial operator- () const;
    /**
     * @brief Substraction of two polynomial.
     * 
     * @param _p the other polynomial
     * @return result
     */
    Polynomial operator- (const Polynomial& _p) const;
    /**
     * @brief Calculate a polynomial minus a real number.
     * 
     * @param _t the real number
     * @return result
     */
    Polynomial operator- (double _t) const;
    /**
     * @brief Calculate a real number minus a polynomial
     * 
     * @param _t the real number
     * @param _p the polynomial 
     * @return result
     */
    friend Polynomial operator- (double _t, const Polynomial& _p);

    /**
     * @brief The multiplication of two polynomial.
     * 
     * @param _p the other polynomial
     * @return result
     */
    Polynomial operator* (const Polynomial& _p) const;
    /**
     * @brief Polynomial multiplies a real number.
     * 
     * @param _t the real number
     * @return result
     */
    Polynomial operator* (double _t) const;
    /**
     * @brief Real number multiplies a polynomial.
     * 
     * @param _t the real number
     * @param _p the polynomial
     * @return result
     */
    friend Polynomial operator* (double _t, const Polynomial& _p);

    /**
     * @brief Get the value of the polynomial at given point.
     * 
     * @param _x given point
     * @return function value
     */
    double operator()(double _x) const;

    /**
     * @brief Output a polynomial.
     * 
     * @param output 
     * @param _p the polynomial
     * @return
     */
    friend std::ostream& operator<< (std::ostream& output, Polynomial& _p);

    /**
     * @brief Input a polynomial.
     * 
     * @param input 
     * @param _p 
     * @return
     */
    friend std::istream& operator>> (std::istream& input, Polynomial& _p);
};
static Polynomial POLY_NULL;

class NewtonInterp
{
private:
    Polynomial interPoly; ///< the interpolation polynomial 
    double2_array tableOfDividedDiffs; ///< table of divided differences

    /**
     * @brief Get the interpolation polynomial from the table of divided differences.
     * 
     * @return 0 for OK
     */
    int getPolyfromTable();

public:
    /**
     * @brief Get the value of the interpolation polynomial at some x.
     * 
     * @param _points interpolation conditions
     * @param _x the given x
     * @return value of interpolation polynomial at _x
     */
    static double Neville_Aitken(const struct InterpConditions _points, double _x);

    /**
     * @brief Create a new table of divided difference according to the given interpolation conditions.
     * 
     * @param _points the interpolation condition
     * @return 0 for OK
     */
    int Newton_overwrite(const struct InterpConditions _points);

    /**
     * @brief Generate a table of divided difference on the base of original table 
     * according to the additional interpolation conditions.
     * 
     * @param _points the additional interpolation conditions
     * @return
     */
    int Newton_addition(const struct InterpConditions _points);

    /**
     * @brief Get the interpolation polynomial.
     * 
     * @return interpolation polynomial.
     */
    Polynomial get_Polynomial();
};

double factorial(int n);

#endif