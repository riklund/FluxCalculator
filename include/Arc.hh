#ifndef Arc_hh
#define Arc_hh 1
#include <iostream>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include "ArcUtils.hh"
#include "VerbosePrinter.hh"
#include "RLException.hh"



using namespace std;

int main(int argc, 
		 char ** argv
		 );

double GetDensityRatio(const vector<double> & xvals, 
					   const char * fileName
					   );

void WriteArcToFile(const vector<double> & xvals, 
					const char * inputFileName, 
					const char * outputFileName,
					double densRatio
					);



#endif
