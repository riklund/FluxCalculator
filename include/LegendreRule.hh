#ifndef LegendreRule_hh
#define LegendreRule_hh 1


#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <ctime>
#include <cstring>
#include <vector>

using namespace std;


/**Computes Gauss-legendre quadrature rules: either the standard one, or a rescaled one.
   Modified from legendre_rule_fast.cpp by John Burkardt, et al, which is availiable under the
   GNU LGPL license. 
   The method used is the so-called Glaser-Liu-Rokhlin method.
 */
class LegendreRule
{
public:
  static vector<pair<double, double> > GetRule(unsigned int numberOfPoints, ///Number of points in the rule.
											   double a = -1, ///The lower bound for rescaling. Default: a=-1, no rescaling.
											   double b = 1 ///Upper bound for rescaling. Default: b=1, no rescaling.
											   ); ///Computes the legendre rule for n points, rescaled to an interval (standard interval is [-1, 1]). Returns a vector with pairs, where the first element of the pair is the position and the second element is the weight. The computations are kind of heavy, so this should ideally be used together with memoization if called frequently.
private:
  LegendreRule() {} ///Thou shall not instantiate this class.


  static void legendre_compute_glr ( int n, double x[], double w[] );
  static void legendre_compute_glr0 ( int n, double *p, double *pp );
  static void legendre_compute_glr1 ( int n, double *roots, double *ders );
  static void legendre_compute_glr2 ( double p, int n, double *roots, double *ders );
  static void legendre_handle ( int n, double a, double b );

  static void rescale ( double a, double b, int n, double x[], double w[] );
  static double rk2_leg ( double t, double tn, double x, int n );
  static double ts_mult ( double *u, double h, int n );

};

#endif
