#ifndef ArcUtils_hh
#define ArcUtils_hh 1

#define EPS 1E-9
#define INIT_VAL (-1337)
#define NUMBER_OF_THETA 5400
#define PI (3.141592653589793238462643)

#include "RLException.hh"
#include <vector>
#include <algorithm>
#include <cmath>
using namespace std;

double InterpolateFlux(const vector<double> & xvals, 
					   const vector<vector<double> > fVals, 
					   double x, 
					   double y
					   );

double TrapezoidIntegrate(const vector<double> & vect, 
						  double delta
						  );

double TrapezoidIntegrate(const vector<vector<double> > & vect,
						  double delta
						  );

double ArcIntegrate(const vector<vector<double> > & vect,
					const vector<double> & xvals
					);
  

#endif
