#include "Density.hh"


////////////
//////////// Code controlling general starting and stopping the program.
////////////

int main(int argc, char *argv[])
{
  CommandLineInterpreter * myInterpreter = InitInterpreter(); 
  try
	{
	  myInterpreter->Initialize(argc, argv);
	  if(!myInterpreter->ReadFlaggedCommand("help").empty() || argc < 1)
		{
		  myInterpreter->PrintHelp();
		  return 0;
		}
	}
  catch(CommandLineException e)
	{
	  cerr << e.what() << endl;
	  return 1;
	}

  ///Read configuration file and act based on that.

  string configFile = myInterpreter->ReadFlaggedCommandStrict("configFile").front().c_str();

  delete myInterpreter;
  
  DensityConfig myConfiguration;
  
  ///If there is no config file, write it and exit. Otherwise, use it.
  if( access(configFile.c_str(), F_OK) == -1 )
	{
	  printf("Config file '%s' not found, aborting.\n", configFile.c_str());
	  return 0;
	}
  else
	{
	  myConfiguration.ReadFile(configFile.c_str());
	}


  ///Initialize stuff.
  
  VerbosePrinter myPrinter(myConfiguration.GetVerbosityLevel());
  myPrinter.Print(3, "Initializing...\n");
  

  PerformCalculations(myConfiguration, myPrinter);
  
  
  myPrinter.Print(2, "Done, exiting.\n");
  return 0;
}


CommandLineInterpreter * InitInterpreter()
{
  CommandLineInterpreter * myInterpreter = new CommandLineInterpreter();

  TMP_LIST_STR(confDefault,"config.conf");
  TMP_LIST_STR(confDescription,"Configuration file");
  myInterpreter->AddCommandLineArgument(CommandLineArgument("configFile", 1, false, "Location of program config file. If the file does not exist, it is created with a default configuration.", confDescription, confDefault));

  myInterpreter->AddCommandLineArgument(CommandLineArgument("help",0,false, "Displays a help message and quits."));

  myInterpreter->SetDescription("Computes probabilities based on output from Compute software.");
  return myInterpreter;
}





////////////
//////////// Calculation code
////////////








void PerformCalculations(const DensityConfig & myConfiguration, VerbosePrinter & myPrinter)
{
  vector<vector<vector<ComplexDouble> > > singleParticleWF(2);
  vector<vector<ComplexDouble> > singleParticleTwiddleFactors(2);
  vector<ComplexDouble> EigendataTwoParticle;
  vector<vector<pair<double, double> > > xGLRules(2);

  double timeStart, timeDelta;
  size_t timeSteps;

  {
	//Basic stuff: curve and weights (needed for renormalization).
	vector<vector<ComplexDouble> > KCurve(2);
	vector<vector<ComplexDouble> > GLWeights(2);
	
	//first element in any eigendata is eigenvalue, rest is eigenvector.
	//for each single particle (outermost vector), we have N eigenvectors (middle vector)
	//containing N elements each. First element in innermost vector is eigenvalue, though.
	vector<vector<vector<ComplexDouble> > > Eigendata(2);
	LoadIndata(myConfiguration, myPrinter, KCurve, GLWeights, Eigendata, EigendataTwoParticle);
	
	InitXGLRules(myConfiguration.GetParticleDomain(), xGLRules);

	timeStart = myConfiguration.GetTimeDomain().start;
	double stop = myConfiguration.GetTimeDomain().stop;
	timeSteps = myConfiguration.GetTimeDomain().precision;
	timeDelta = (stop - timeStart) / timeSteps;

	myPrinter.Print(3, "Initializing twiddle factors...");
	FillSingleParticleTwiddle(Eigendata, timeDelta, myConfiguration.GetUnits().hbarTimesLambda, singleParticleTwiddleFactors);
	myPrinter.Print(3, "done.\n");
	myPrinter.Print(3, "Computing single-particle spatial wavefunctions...");	
	FillSingleParticleSpatial(xGLRules, KCurve, GLWeights, Eigendata, singleParticleWF);
	myPrinter.Print(3, "done.\n");

  }

  
  myPrinter.Print(3, "Opening output files.\n");	
  ManagedOutputFiles myFiles(myPrinter, myConfiguration.GetOutputFiles());
  vector<vector<ComplexDouble> > twoParticleWF(xGLRules.at(0).size(), vector<ComplexDouble>(xGLRules.at(1).size()));
  myPrinter.Print(3, "Done opening output files.\n");
  
  myPrinter.Print(3, "Starting time stepping with stepsize %5.2e.\n", timeDelta);

  double twoParticleNorm = -1337;


  for(size_t timeIndex = 0; timeIndex < timeSteps; ++timeIndex)
	{
	  double time = timeStart + timeDelta * timeIndex;
	  myPrinter.Print(2, "Time is now %+5.2e, ( %d / %d steps ).\n", time, timeIndex, timeSteps);
	  myPrinter.Print(5, "Saving single particle densities...");
	  for(uint i = 0; i<2; ++i)
		{
		  size_t bState = myConfiguration.GetOutputFiles().SingleBasis.at(i);
		  double oneRho = ComputeSingleParticleRho(xGLRules.at(i), singleParticleWF.at(i).at(bState));
		  SaveRhoToFile(myFiles.SingleBasis.at(i), time, oneRho);
		}
	  myPrinter.Print(5, "done.\n");
	  myPrinter.Print(5, "Creating two-particle wavefunction...");
	  FillTwoParticleWF(singleParticleWF, EigendataTwoParticle, twoParticleWF);
	  myPrinter.Print(5, "done.\n");
	  myPrinter.Print(5, "Computing two-particle density and saving...");
	  double twoRho = ComputeTwoParticleRho(xGLRules, twoParticleWF);
	  if(twoParticleNorm > -100)
		{
		  twoRho /= twoParticleNorm;
		}
	  else
		{
		  twoParticleNorm = twoRho;
		  twoRho = 1;
		}
	  SaveRhoToFile(myFiles.TimeSeries, time, twoRho);
	  myPrinter.Print(5, "done.\n");

	  myPrinter.Print(5, "Evolving time...");
	  EvolveTime(singleParticleTwiddleFactors, singleParticleWF);
	  myPrinter.Print(5, "done.\n");
	}
}

void EvolveTime(const vector<vector<ComplexDouble> > & singleParticleTwiddleFactors, vector<vector<vector<ComplexDouble> > > & singleParticleWF)
{
  ///Evolve time.
  for(uint i = 0; i<2; ++i)
	{
#pragma omp parallel for
	  for(size_t j = 0; j<singleParticleWF.at(i).size(); ++j)
		{
		  for(size_t x = 0; x<singleParticleWF.at(i).at(j).size(); ++x)
			{
			  singleParticleWF.at(i).at(j).at(x) *= singleParticleTwiddleFactors.at(i).at(j);
			}
		}
	}  
}


void InitXGLRules(const vector<DomainSpecific> & particleDomain, vector<vector<pair<double, double> > > & xGLRules)
{
  for(uint i = 0; i<2; ++i)
	{
	  double start = particleDomain.at(i).start;
	  double stop = particleDomain.at(i).stop;
	  size_t precision = particleDomain.at(i).precision;
	  xGLRules.at(i) = LegendreRule::GetRule(precision, start, stop);
	}
}


double ComputeSingleParticleRho(const vector<pair<double, double> > & xGLRule, const vector<ComplexDouble> & singleParticleWF)
{
  if(xGLRule.size() != singleParticleWF.size())
	throw RLException("Size mismatch: xGLRule had size %d, singleParticleWF had size %d.", xGLRule.size(), singleParticleWF.size());

  double density = 0;
  for(size_t x = 0; x<xGLRule.size(); ++x)
	{
	  density += xGLRule.at(x).second * pow(abs(singleParticleWF.at(x)),2.0);
	}
  return density;
}


void LoadIndata(const DensityConfig & myConfiguration, VerbosePrinter & myPrinter, vector<vector<ComplexDouble> > & KCurve, vector<vector<ComplexDouble> > & GLWeights, vector<vector<vector<ComplexDouble> > > &Eigendata, vector<ComplexDouble> & EigendataTwoParticle)
{
  myPrinter.Print(2, "Loading input data.\n");
  for(uint i = 0; i<2; ++i)
	{
	  LoadFileToVector(myConfiguration.GetInputFiles().KCurveParticle.at(i).c_str(), myPrinter, KCurve[i]);
	  LoadFileToVector(myConfiguration.GetInputFiles().GLWeightsParticle.at(i).c_str(), myPrinter, GLWeights[i]);
	  LoadEigenInfo(myConfiguration.GetInputFiles().EigendataParticle.at(i).c_str(), myPrinter, Eigendata[i]);
	}
  vector<vector<ComplexDouble> > tempEigenTwo;
  LoadEigenInfo(myConfiguration.GetInputFiles(). EigendataTwoParticles.c_str(), myPrinter, tempEigenTwo);
  if(tempEigenTwo.size() != 1)
	{
	  throw RLException("Found %d two-particle eigenvectors, one expected.", tempEigenTwo.size());
	}
  EigendataTwoParticle = tempEigenTwo.at(0);


  //Assert stuff.
  if(KCurve[0].size() != KCurve[1].size())
	throw RLException("Unexpected difference in single-particle basis size: KCurve ");
  if(GLWeights[0].size() != GLWeights[1].size())
	throw RLException("Unexpected difference in single-particle basis size: GLWeights (%d vs %d).", GLWeights[0].size(), GLWeights[1].size());
  if(Eigendata[0].size() != Eigendata[1].size())
	throw RLException("Unexpected difference in single-particle basis size: Eigendata");
  if(KCurve[0].size() * 2 != Eigendata[0].size())
	throw RLException("Unexpected difference between eigendata (%d) and KCurve (%d) size.", Eigendata[0].size(), KCurve[0].size());
  if(KCurve[0].size() != GLWeights[0].size())
	throw RLException("Unexpected difference between KCurve and GLWeights size.");
  if((Eigendata[0].size())*(Eigendata[1].size()) != EigendataTwoParticle.size()-1)
	throw RLException("Size mismatch: single particle %d, two-particle %d", Eigendata[0].size(), EigendataTwoParticle.size());

  myPrinter.Print(2, "Done loading data.\n");
}

void LoadFileToVector(string fileName, VerbosePrinter & myPrinter, vector<ComplexDouble> & output)
{
  myPrinter.Print(3,"Loading file '%s'...", fileName.c_str());

  output.clear();
  ifstream fin(fileName.c_str(), fstream::in);
  double re, im;
  while(fin >> re >> im)
	{
	  output.push_back(ComplexDouble(re, im));
	}

  myPrinter.Print(3,"done.\n");
}


void LoadEigenInfo(string fileName, VerbosePrinter & myPrinter, vector<vector<ComplexDouble> > & output)
{
  myPrinter.Print(3,"Loading file '%s'...", fileName.c_str());
  output.clear();
  
  ifstream fin(fileName.c_str(), fstream::in);
  if(!fin.good())
	{
	  throw RLException("Failure on opening file '%s'.", fileName.c_str());
	}

  string current;

  while(!(getline(fin, current).eof()))
	{
	  if(current.find("#")!=string::npos)
		continue;
	  output.push_back(vector<ComplexDouble>());
	  stringstream mystream(current);
	  double re, im;
	  while(mystream >> re >> im)
		output.back().push_back(ComplexDouble(re, im));
	}

  for(size_t i = 1; i<output.size(); ++i)
	{
	  if(output.at(i).size() != output.at(i-1).size())
		{
		  throw RLException("Read dimension mismatch.");
		}
	}
  
  myPrinter.Print(3,"done.\n");
}


void FillTwoParticleWF(const vector<vector<vector<ComplexDouble > > > & singleParticleWF, const vector<ComplexDouble> & EigendataTwoParticle, vector<vector<ComplexDouble> > & twoParticleWF)

{
#pragma omp parallel for
  for(size_t xa = 0; xa<twoParticleWF.size(); ++xa)
	{
	  for(size_t xb = 0; xb<twoParticleWF.at(xa).size(); ++xb)
		{
		  twoParticleWF.at(xa).at(xb) = TwoParticleWF(xa, xb, singleParticleWF, EigendataTwoParticle);
		}
	}
}

ComplexDouble TwoParticleWF(size_t xa, size_t xb, const vector<vector<vector<ComplexDouble> > > & singleParticleWF, const vector<ComplexDouble> & EigendataTwoParticle)
{
  ComplexDouble toReturn = 0.0;
  size_t N1 = singleParticleWF.at(0).size();
  size_t N2 = singleParticleWF.at(1).size();
  size_t Ntot = N1 * N2;
  for(size_t i = 0; i<Ntot; ++i)
	{
	  size_t a = i / N1;
	  size_t b = i % N2;
	  toReturn += EigendataTwoParticle.at(i+1) * singleParticleWF.at(0).at(a).at(xa) * singleParticleWF.at(1).at(b).at(xb);
	}
  return toReturn;
}



void FillSingleParticleSpatial(const vector<vector<pair<double, double> > > & xGLRules, const vector<vector<ComplexDouble> > & KCurve, const vector<vector<ComplexDouble> > & GLWeights, const vector<vector<vector<ComplexDouble> > > & Eigendata, vector<vector<vector<ComplexDouble> > > & singleParticleSpatialWF)
{
  singleParticleSpatialWF.resize(2);
  for(uint i = 0; i<2; ++i)
	{
	  singleParticleSpatialWF.at(i).resize(Eigendata.at(i).size(), vector<ComplexDouble>(xGLRules.at(i).size(), 0.0));
	}
  
  for(uint partIndex = 0; partIndex < 2; ++partIndex)
	{
#pragma omp parallel for
	  for(size_t x = 0; x<xGLRules.at(partIndex).size(); ++x)
		{
		  for(size_t eigIndex = 0; eigIndex < Eigendata.at(partIndex).size(); ++eigIndex)
			{
			  singleParticleSpatialWF.at(partIndex).at(eigIndex).at(x) = 
				SingleParticleSpatial(KCurve.at(partIndex), 
									  GLWeights.at(partIndex), 
									  Eigendata.at(partIndex).at(eigIndex), 
									  xGLRules.at(partIndex).at(x).first
									  );
			}
		}
	}
}


ComplexDouble SingleParticleSpatial(const vector<ComplexDouble> & KCurve, const vector<ComplexDouble> GLWeights, const vector<ComplexDouble> & Eigendata, double x)
{
  ComplexDouble toReturn = 0.0;

  vector<ComplexDouble (*)(const ComplexDouble&)> bfun(2);
  bfun.at(0) = &std::sin, bfun.at(1) = &std::cos;

  size_t N = KCurve.size();

  for(size_t i = 0; i<2*N; ++i)
	{

	  size_t curvePointer = i % N;
	  size_t basisPointer = i / N;

	  toReturn += 
		//coefficient
		Eigendata.at(i+1) *
		//spatial part
		sqrt(GLWeights.at(curvePointer)) *
		sqrt(1/PI) * 
		bfun.at(basisPointer)(KCurve.at(curvePointer) * x);
	}

  return toReturn;
}


void FillSingleParticleTwiddle(const vector<vector<vector<ComplexDouble> > > & Eigendata, const double timestep, const double hbar, vector<vector<ComplexDouble> > & singleParticleTwiddleFactors)
{
  for(uint partIndex = 0; partIndex<2; ++partIndex)
	{
	  singleParticleTwiddleFactors.at(partIndex).resize(Eigendata.at(partIndex).size(), 0.0);
	  for(size_t i = 0; i<Eigendata.at(partIndex).size(); ++i)
		{
		  singleParticleTwiddleFactors.at(partIndex).at(i) = 
			PsiSingleTwiddle(Eigendata.at(partIndex).at(i).at(0), timestep, hbar);
		}
	}
}


ComplexDouble PsiSingleTwiddle(ComplexDouble eigenvalue, double timestep, double hbar)
{
  ComplexDouble toReturn = exp(ComplexDouble(0,-1.0)*
				  eigenvalue *
				  timestep / hbar);
  return toReturn;
}

double ComputeTwoParticleRho(const vector<vector<pair<double, double> > > & xGLRules, const vector<vector<ComplexDouble> > twoParticleWF)
{
  double rho = 0;
  for(size_t xa = 0; xa < xGLRules.at(0).size(); ++xa)
	{
	  for(size_t xb = 0; xb < xGLRules.at(1).size(); ++xb)
		{
		  rho += 
			xGLRules.at(0).at(xa).second * 
			xGLRules.at(1).at(xb).second *
			pow(abs(twoParticleWF.at(xa).at(xb)), 2);
		}
	}
  return rho;
}

void SaveRhoToFile(FILE * fout, double time, double rho)
{
  if(fout)
	{
	  fprintf(fout, "%+13.10e %+13.10e\n", time, rho);
	  fflush(fout);
	}
}
