#ifndef ComputeConfig_hh
#define ComputeConfig_hh 1

#include <list>
#include <vector>

#include "libconfig.h++"
#include "RLException.hh"
#include "RLMacros.hpp"

using namespace std;
using namespace libconfig;

#define CONFIG_FILE_VERSION 0.10

class InputFiles
{
public:
  InputFiles()
	:KCurveParticle(2, ""), GLWeightsParticle(2, ""), EigendataParticle(2, ""), EigendataTwoParticles("")
  {}
  vector<string> KCurveParticle;
  vector<string> GLWeightsParticle;
  vector<string> EigendataParticle;
  
  string EigendataTwoParticles;
};

class DomainSpecific
{
public:
  DomainSpecific()
	:start(0.0), stop(0.0), precision(0)
  { }
  double start;
  double stop;
  size_t precision;
};

class UnitSet
{
public:
  UnitSet()
	:hbarTimesLambda(1), massOverLambda2(1)
  { }
  double hbarTimesLambda;
  double massOverLambda2;
};

class FluxConfig
{
public:
  FluxConfig();
  FluxConfig(const char * fileName
				);
  ~FluxConfig();
  void ReadFile(const char * fileName
				);


  //Get/set
  uint GetVerbosityLevel() const;
  uint GetNumberOfThreads() const;
  const InputFiles & GetInputFiles() const;
  const vector<DomainSpecific> & GetParticleDomain() const;
  const UnitSet & GetUnits() const;
  const vector<size_t> & GetSingleParticleBasisNumber() const;

  const string & GetOutputDensityFile() const;
  const string & GetOutputFluxFile() const;

protected:
  void InitProgramGenerals(Setting & root);
  void InitInputFiles(Setting & root);
  void InitDomainSpecifics(Setting & root);
  void InitUnits(Setting & root);
  void InitSingleParticleBasisNumber(Setting & root);

  static void SizetLookup(const Setting & root, const char * propName, size_t & output);

private:
  uint verbosityLevel;
  uint numberOfThreads;
  
  InputFiles myInputFiles;

  vector<size_t> mySingleParticleBasisNumber;
  vector<DomainSpecific> particleDomain;

  string outputDensity;
  string outputFlux;

  UnitSet myUnits;
};

#endif
