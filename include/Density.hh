#ifndef Density_hh
#define Density_hh 1

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <unistd.h>
#include <algorithm>
#include <fstream>
#include <functional>

#include "DensityConfig.hh"
#include "CommandLineInterpreter.hh"
#include "VerbosePrinter.hh"
#include "RLMacros.hpp"
#include "LegendreRule.hh"
#include "ManagedOutputFiles.hh"


using namespace std;

//Used to init command line arguments.
#define TMP_LIST_STR(n, s) list<string> n; n.push_back(s);



int main(int argc, 
		 char *argv[]
		 );

CommandLineInterpreter * InitInterpreter();///Returns a command line interpreter with defined commands.


void PerformCalculations(const DensityConfig & myConfiguration, 
						 VerbosePrinter & myPrinter
						 );

void InitXGLRules(const vector<DomainSpecific> & particleDomain, 
				  vector<vector<pair<double, double> > > & xGLRules
				  );

double ComputeSingleParticleRho(const vector<pair<double, double> > & xGLRule, 
								const vector<ComplexDouble> & singleParticleWF
								);

void LoadIndata(const DensityConfig & myConfiguration, 
				VerbosePrinter & myPrinter, 
				vector<vector<ComplexDouble> > & KCurve, 
				vector<vector<ComplexDouble> > & GLWeights, 
				vector<vector<vector<ComplexDouble> > > &Eigendata, 
				vector<ComplexDouble> & EigendataTwoParticle
				);

void LoadFileToVector(string fileName, 
					  VerbosePrinter & myPrinter, 
					  vector<ComplexDouble> & output
					  );

void LoadEigenInfo(string fileName, 
				   VerbosePrinter & myPrinter, 
				   vector<vector<ComplexDouble> > & output
				   );

void FillTwoParticleWF(const vector<vector<vector<ComplexDouble > > > & singleParticleWF, 
					   const vector<ComplexDouble> & EigendataTwoParticle,
					   vector<vector<ComplexDouble> > & twoParticleWF
					   );

ComplexDouble TwoParticleWF(size_t xa, 
							size_t xb, 
							const vector<vector<vector<ComplexDouble> > > & singleParticleWF,
							const vector<ComplexDouble> & EigendataTwoParticle
							);

void FillSingleParticleSpatial(const vector<vector<pair<double, double> > > & xGLRules, 
							   const vector<vector<ComplexDouble> > & KCurve, 
							   const vector<vector<ComplexDouble> > & GLWeights, 
							   const vector<vector<vector<ComplexDouble> > > & Eigendata, 
							   vector<vector<vector<ComplexDouble> > > & singleParticleSpatialWF
							   );


ComplexDouble SingleParticleSpatial(const vector<ComplexDouble> & KCurve, 
									const vector<ComplexDouble> GLWeights, 
									const vector<ComplexDouble> & Eigendata, 
									const double x
									);

void FillSingleParticleTwiddle(const vector<vector<vector<ComplexDouble> > > & Eigendata, 
							   const double timestep, 
							   const double hbar,
							   vector<vector<ComplexDouble> > & singleParticleTwiddleFactors
							   );

ComplexDouble PsiSingleTwiddle(ComplexDouble eigenvalue, 
							   double timestep, 
							   double hbar
							   );

double ComputeTwoParticleRho(const vector<vector<pair<double, double> > > & xGLRules, 
						   const vector<vector<ComplexDouble> > twoParticleWF
						   );

void SaveRhoToFile(FILE * fout, 
				   double time, 
				   double rho
				   );


void EvolveTime(const vector<vector<ComplexDouble> > & singleParticleTwiddleFactors, 
				vector<vector<vector<ComplexDouble> > > & singleParticleWF
				);

#endif
