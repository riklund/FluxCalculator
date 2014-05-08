#ifndef Compute_hh
#define Compute_hh 1

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <unistd.h>
#include <algorithm>

#include "DensityConfig.hh"
#include "CommandLineInterpreter.hh"
#include "VerbosePrinter.hh"
#include "RLMacros.hpp"
#include "LegendreRule.hh"

using namespace std;

//Used to init command line arguments.
#define TMP_LIST_STR(n, s) list<string> n; n.push_back(s);


int main(int argc, 
		 char *argv[]
		 );


CommandLineInterpreter * InitInterpreter();///Returns a command line interpreter with defined commands.


void PerformCalculations(const DensityConfig & myConfiguration, VerbosePrinter & myPrinter);

ComplexDouble PsiSingle(uint partIndex, size_t x, size_t t, size_t eigIndex);

ComplexDouble PsiTwo(size_t xa, size_t xb, size_t t);

void LoadIndata(const DensityConfig & myConfiguration, VerbosePrinter & myPrinter);

void LoadFileToVector(string fileName, VerbosePrinter & myPrinter, vector<ComplexDouble> & output);

void FillUnits(const DensityConfig & myConfiguration, VerbosePrinter & myPrinter);

void FillPsiTwoCache(const DensityConfig & myConfiguration, VerbosePrinter & myPrinter);

void Fill3DGrid(const DensityConfig & myConfiguration, VerbosePrinter & myPrinter);

void FillPsiSingleCache(const DensityConfig & myConfiguration, VerbosePrinter & myPrinter);

ComplexDouble PsiSingle(uint partIndex, size_t x, size_t t, size_t eigIndex);

ComplexDouble KToE(ComplexDouble k);

ComplexDouble PsiTwo(size_t xa, size_t xb, size_t t);

void FillRhoT(const DensityConfig & myConfiguration, VerbosePrinter & myPrinter);

void SaveRhoTToFile(const DensityConfig & myConfiguration, VerbosePrinter & myPrinter);

///Global variables: used for calculations.

//This is our 3D grid.
//We need to integrate over the x:es
vector<vector<pair<double, double> > > xGLRules(2);
//but not over the t:s
vector<double> tValues;


//Basic stuff: curve and weights (needed for renormalization).
vector<vector<ComplexDouble> > KCurve(2);
vector<vector<ComplexDouble> > GLWeights(2);

//first element in any eigendata is eigenvalue, rest is eigenvector.
//for each single particle (outermost vector), we have N eigenvectors (middle vector)
//containing N elements each. First element in innermost vector is eigenvalue, though.
vector<vector<ComplexDouble> > Eigendata(2);

vector<double> rhoT;
double hbarTimesLambda = 0, massOverLambda2 = 0;

//For two particles, we only have one eigenvector: the one of the resonance.
//First element is eigenvalue, rest is eigenvector.
vector<ComplexDouble> EigendataTwoParticle;

///From outer to inner:
/// Particle ID
/// arg 1
/// arg 2 
/// arg 3
vector<vector<vector<vector<ComplexDouble> > > > psiSingleCache;

///From outer to inner:
///arg 1, arg 2, arg 3
vector<vector<vector<ComplexDouble> > > psiTwoCache;

#endif
