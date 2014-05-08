#include "Density.hh"

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


void PerformCalculations(const DensityConfig & myConfiguration, VerbosePrinter & myPrinter)
{  
  LoadIndata(myConfiguration, myPrinter);
  FillUnits(myConfiguration, myPrinter);
  Fill3DGrid(myConfiguration, myPrinter);
  FillPsiSingleCache(myConfiguration, myPrinter);
  FillPsiTwoCache(myConfiguration, myPrinter);
  FillRhoT(myConfiguration, myPrinter);
  SaveRhoTToFile(myConfiguration, myPrinter);
}

void LoadIndata(const DensityConfig & myConfiguration, VerbosePrinter & myPrinter)
{
  myPrinter.Print(2, "Loading input data.\n");
  for(uint i = 0; i<2; ++i)
	{
	  LoadFileToVector(myConfiguration.GetInputFiles().KCurveParticle.at(i).c_str(), myPrinter, KCurve[i]);
	  LoadFileToVector(myConfiguration.GetInputFiles().GLWeightsParticle.at(i).c_str(), myPrinter, GLWeights[i]);
	  LoadFileToVector(myConfiguration.GetInputFiles().EigendataParticle.at(i).c_str(), myPrinter, Eigendata[i]);
	}

  LoadFileToVector(myConfiguration.GetInputFiles(). EigendataTwoParticles.c_str(), myPrinter, EigendataTwoParticle);



  //Assert stuff.
  if(KCurve[0].size() != KCurve[1].size())
	throw RLException("Unexpected difference in single-particle basis size: KCurve ");
  if(GLWeights[0].size() != GLWeights[1].size())
	throw RLException("Unexpected difference in single-particle basis size: GLWeights ");
  if(Eigendata[0].size() != Eigendata[1].size())
	throw RLException("Unexpected difference in single-particle basis size: Eigendata");
  if(KCurve[0].size() + 1 != Eigendata[0].size())
	throw RLException("Unexpected difference between eigendata and KCurve size.");
  if(KCurve[0].size() != GLWeights[0].size())
	throw RLException("Unexpected difference between KCurve and GLWeights size.");

}


void FillUnits(const DensityConfig & myConfiguration, VerbosePrinter & myPrinter)
{
  hbarTimesLambda = myConfiguration.GetUnits().hbarTimesLambda;
  massOverLambda2 = myConfiguration.GetUnits().massOverLambda2;
}

void LoadFileToVector(string fileName, VerbosePrinter & myPrinter, vector<ComplexDouble> & output)
{
  myPrinter.Print(3,"Loading file %s...", fileName.c_str());
  FILE * infile = fopen(fileName.c_str(), "r");
  if(infile == NULL)
	{
	  throw RLException("Could not open input file '%s'.",fileName.c_str());
	}
  double re, im;
  while(fscanf(infile, "%lf %lf\n", &re, &im) != EOF)
	{
	  output.push_back(ComplexDouble(re, im));
	}
  fclose(infile);
  myPrinter.Print(3,"done.\n");
}

void FillPsiTwoCache(const DensityConfig & myConfiguration, VerbosePrinter & myPrinter)
{
  psiTwoCache.resize(xGLRules.at(0).size(), vector<vector<ComplexDouble> >(xGLRules.at(1).size(), vector<ComplexDouble>(tValues.size(), 0.0)));
  for(size_t xa = 0; xa<xGLRules.at(0).size(); ++xa)
	{
#pragma omp parallel for
	  for(size_t xb = 0; xb<xGLRules.at(1).size(); ++xb)
		{
		  for(size_t t = 0; t<tValues.size(); ++t)
			{
			  psiTwoCache.at(xa).at(xb).at(t) = PsiTwo(xa, xb, t);
			}
		}
	}
}

void Fill3DGrid(const DensityConfig & myConfiguration, VerbosePrinter & myPrinter)
{
  for(uint i = 0; i<2; ++i)
	{
	  double start = myConfiguration.GetParticleDomain().at(i).start;
	  double stop = myConfiguration.GetParticleDomain().at(i).stop;
	  size_t precision = myConfiguration.GetParticleDomain().at(i).precision;

	  ///Construct GL rule.
	  for(size_t j = 0; j<precision; ++j)
		{
		  xGLRules.at(i) = LegendreRule::GetRule(precision, start, stop);
		}
	}

  double start = myConfiguration.GetTimeDomain().start;
  double stop = myConfiguration.GetTimeDomain().stop;
  size_t precision = myConfiguration.GetTimeDomain().precision;
  tValues.resize(precision);
  double delta = (stop - start) / precision;
  for(size_t i = 0; i<precision; ++i)
	{
	  tValues.at(i) = start + delta * i;
	}
}
						
void FillPsiSingleCache(const DensityConfig & myConfiguration, VerbosePrinter & myPrinter)
{
  psiSingleCache.resize(2);
  for(uint i = 0; i<2; ++i)
	psiSingleCache.at(i).resize(xGLRules.at(i).size(), vector<vector<ComplexDouble> >(tValues.size(), vector<ComplexDouble>(Eigendata.at(i).size(), 0) ) );
  
  for(uint partIndex = 0; partIndex < 2; ++partIndex)
	{
#pragma omp parallel for
	  for(size_t x = 0; x<xGLRules.at(partIndex).size(); ++x)
		{
		  for(size_t t = 0; t<tValues.size(); ++t)
			{
			  for(size_t eigIndex = 0; eigIndex < Eigendata.at(partIndex).size(); ++eigIndex)
				{
				  psiSingleCache.at(partIndex).at(x).at(t).at(eigIndex) = 
					PsiSingle(partIndex, x, t, eigIndex);
				}
			}
		}
	}
}


ComplexDouble PsiSingle(uint partIndex, size_t x, size_t t, size_t eigIndex)
{
  ComplexDouble toReturn = 0.0;

  vector<ComplexDouble (*)(const ComplexDouble&)> bfun(2);
  bfun.at(0) = &std::sin, bfun.at(1) = &std::cos;

  ///sine functions.
  for(size_t i = 0; i<KCurve.at(partIndex).size(); ++i)
	{
	  double XX = xGLRules.at(partIndex).at(x).first;

	  toReturn += 
		Eigendata.at(partIndex).at(i+1) *
		sqrt(GLWeights.at(partIndex).at(i)) *
		sqrt(1/PI) * 
		bfun.at(i%2)(KCurve.at(partIndex).at(i) * XX) *
		exp(ComplexDouble(0,-1.0)*KToE(KCurve.at(partIndex).at(i)) * tValues.at(t) / hbarTimesLambda);
	}
  return toReturn;
}

ComplexDouble KToE(ComplexDouble k)
{
  return pow(hbarTimesLambda*k,2)/(2*massOverLambda2);
}
 
ComplexDouble PsiTwo(size_t xa, size_t xb, size_t t)
{
  ComplexDouble toReturn = 0.0;
  size_t N1 = KCurve.at(0).size();
  size_t N2 = KCurve.at(1).size();
  size_t Ntot = N1 * N2;
  for(size_t i = 0; i<Ntot; ++i)
	{
	  size_t a = i / N1;
	  size_t b = i % N2;
	  toReturn += psiSingleCache[0][xa][t][a] * psiSingleCache[1][xb][t][b];
	}
  return toReturn;
}



void FillRhoT(const DensityConfig & myConfiguration, VerbosePrinter & myPrinter)
{
  rhoT.resize(tValues.size(), 0.0);
  for(size_t t = 0; t<tValues.size(); ++t)
	{
	  for(size_t xa = 0; xa < xGLRules.at(0).size(); ++xa)
		{
		  for(size_t xb = 0; xb < xGLRules.at(1).size(); ++xb)
			{
			  rhoT.at(t) += 
				xGLRules.at(0).at(xa).second * 
				xGLRules.at(1).at(xb).second *
				pow(abs(psiTwoCache.at(xa).at(xb).at(t)), 2);
			}
		}
	}
}



void SaveRhoTToFile(const DensityConfig & myConfiguration, VerbosePrinter & myPrinter)
{
  string fileName = myConfiguration.GetOutputFiles().TimeSeries;
  FILE * fout = fopen(fileName.c_str(), "w");
  if(fout == NULL)
	{
	  throw RLException("Could not open output file '%s'.", fileName.c_str());
	}
  for(size_t i = 0; i < rhoT.size(); ++i)
	{
	  fprintf(fout, "%+13.10e %+13.10e\n", tValues.at(i), rhoT.at(i));
	}
  
  fclose(fout);
}
