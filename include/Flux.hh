#ifndef Flux_hh
#define Flux_hh 1

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

#include "FluxConfig.hh"
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


void SaveTwoParticleDensity(const vector<vector<pair<double, double> > > & xGLRules, const   vector<vector<ComplexDouble> > & twoParticleWF, const string & filename);


void PerformCalculations(const FluxConfig & myConfiguration, 
						 VerbosePrinter & myPrinter
						 );

void InitXGLRules(const vector<DomainSpecific> & particleDomain, 
				  vector<vector<pair<double, double> > > & xGLRules,
				  bool useGLRule
				  );

vector<pair<double, double> > GetNormalRule(size_t precision, double start, double stop);

double ComputeSingleParticleRho(const vector<pair<double, double> > & xGLRule, 
								const vector<ComplexDouble> & singleParticleWF
								);

void LoadIndata(const FluxConfig & myConfiguration, 
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

void FillTwoParticleWFGrad(const vector<vector<vector<ComplexDouble > > > & singleParticleWF, 
						   const vector<vector<vector<ComplexDouble> > >& singleParticleWFPrime, 
						   const vector<ComplexDouble> & EigendataTwoParticle,
						   vector<vector<pair<ComplexDouble, ComplexDouble> > > & twoParticleWFGrad
					   );

ComplexDouble TwoParticleWF(size_t xa, 
							size_t xb, 
							const vector<vector<vector<ComplexDouble> > > & singleParticleWF,
							const vector<ComplexDouble> & EigendataTwoParticle
							);

pair<ComplexDouble, ComplexDouble> TwoParticleWFGrad(size_t xa, 
								 size_t xb, 
								 const vector<vector<vector<ComplexDouble> > > & singleParticleWF,
								 const vector<vector<vector<ComplexDouble> > > & singleParticleWFPrime,
								 const vector<ComplexDouble> & EigendataTwoParticle
								 );

void FillSingleParticleSpatial(const vector<vector<pair<double, double> > > & xGLRules, 
							   const vector<vector<ComplexDouble> > & KCurve, 
							   const vector<vector<ComplexDouble> > & GLWeights, 
							   const vector<vector<vector<ComplexDouble> > > & Eigendata, 
							   vector<vector<vector<ComplexDouble> > > & singleParticleSpatialWF, 
							   const bool diff = false
							   );


ComplexDouble SingleParticleSpatial(const vector<ComplexDouble> & KCurve, 
									const vector<ComplexDouble> GLWeights, 
									const vector<ComplexDouble> & Eigendata, 
									const double x,
									const bool diff = false
									);


double ComputeTwoParticleRho(const vector<vector<pair<double, double> > > & xGLRules, 
						   const vector<vector<ComplexDouble> > twoParticleWF
						   );

void CoutEqualBars(size_t count);

pair<ComplexDouble, ComplexDouble> operator*(const pair<ComplexDouble, ComplexDouble> & p, ComplexDouble val);

pair<ComplexDouble, ComplexDouble> operator-(const pair<ComplexDouble, ComplexDouble> & left, const pair<ComplexDouble, ComplexDouble> & right);

pair<ComplexDouble, ComplexDouble> conj(const pair<ComplexDouble, ComplexDouble> & p);

#endif
