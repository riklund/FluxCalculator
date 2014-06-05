#include "FluxConfig.hh"

FluxConfig::FluxConfig()
  :verbosityLevel(0), numberOfThreads(1), particleDomain(2), outputDensity("")
{
  
}

FluxConfig::FluxConfig(const char * fileName)
 : particleDomain(2)
{
  ReadFile(fileName);
}

void FluxConfig::ReadFile(const char * fileName)
{
  Config cfg;
  try
	{
	  cfg.readFile(fileName);
	}
  catch(ParseException &ex)
	{
	  printf("Caught ParseException : %s Location : line %d in the config file.\n", ex.getError(), ex.getLine());
	  throw ex;
	}
  Setting & root = cfg.getRoot();
  

  double version;
  if( root.lookupValue("Version", version) )
	{
	  if( ! DBL_EQUAL(version, CONFIG_FILE_VERSION) )
		throw RLException("Invalid configuration file version: was %f, required %f", version, CONFIG_FILE_VERSION);
	}
  else
	{
	  throw RLException("Could not find property 'Version' in configuration file.");
	}
  
  InitProgramGenerals(root);
  InitSingleParticleBasisNumber(root);
  InitInputFiles(root);
  InitDomainSpecifics(root);
  InitUnits(root);
}

void FluxConfig::InitProgramGenerals(Setting & root)
{
  if( ! root.exists("Program") || !root["Program"].isGroup())
	throw RLException("Could not look up group 'Program' in config file.");
  Setting & program = root["Program"];
  if( ! program.lookupValue("VerbosityLevel", verbosityLevel) )
	throw RLException("Could not look up verbosity level.");
  if( ! program.lookupValue("NumberOfThreads", numberOfThreads) )
	throw RLException("Could not look up number of threads.");
}

void FluxConfig::InitInputFiles(Setting & root)
{
  if( ! root.exists("InputFiles") || ! root["InputFiles"].isGroup() )
	throw RLException("Invalid input file group");
  Setting & input = root["InputFiles"];
  if( ! input.lookupValue("KCurveParticleA", myInputFiles.KCurveParticle.at(0)) )
	throw RLException("Could not look up value for KCurveParticleA input file.");
  if( ! input.lookupValue("KCurveParticleB", myInputFiles.KCurveParticle.at(1)) )
	throw RLException("Could not look up value for KCurveParticleB input file.");
  if( ! input.lookupValue("GLWeightsParticleA", myInputFiles.GLWeightsParticle.at(0)) )
	throw RLException("Could not look up value for GLWeightsParticleA input file.");
  if( ! input.lookupValue("GLWeightsParticleB", myInputFiles.GLWeightsParticle.at(1)) )
	throw RLException("Could not look up value for GLWeightsParticleB input file.");

  if( ! input.lookupValue("EigendataParticleA", myInputFiles.EigendataParticle.at(0)) )
	throw RLException("Could not look up value for EigendataParticleA input file.");
  if( ! input.lookupValue("EigendataParticleB", myInputFiles.EigendataParticle.at(1)) )
	throw RLException("Could not look up value for EigendataParticleB input file.");
  if( ! input.lookupValue("EigendataTwoParticles", myInputFiles.EigendataTwoParticles) )
	throw RLException("Could not look up value for EigendataTwoParticles input file.");
}

void FluxConfig::InitSingleParticleBasisNumber(Setting & root)
{
 if( ! root.exists("Output") || ! root["Output"].isGroup() )
	throw RLException("Invalid SingleBasis file group");
 Setting & basis = root["Output"];
 mySingleParticleBasisNumber.resize(2, 0);
 int temp;
 if( ! basis.lookupValue("SingleA", temp))
   throw RLException("Could not look up SingleA value.");
 mySingleParticleBasisNumber.at(0) = temp;
 if( ! basis.lookupValue("SingleB", temp))
   throw RLException("Could not look up SingleB value.");
 mySingleParticleBasisNumber.at(1) = temp;
 if( ! basis.lookupValue("DensityFile", outputDensity))
   throw RLException("Could not look up output density.");
 if( ! basis.lookupValue("FluxFile", outputFlux))
   throw RLException("Could not look up output flux.");
}


void FluxConfig::SizetLookup(const Setting & root, const char * propName, size_t & output)
{
  unsigned int tmp;
  if(!root.lookupValue(propName, tmp))
	{
	  throw RLException("Lookup of property '%s' failed.", propName);
	}
  output=tmp;
  return;
}

void FluxConfig::InitDomainSpecifics(Setting & root)
{
  if( ! root.exists("DomainSpecifics") || ! root["DomainSpecifics"].isGroup() )
	throw RLException("Invalid DomainSpecifics file group");
  Setting & dspec = root["DomainSpecifics"];

  if( ! dspec.exists("ParticleANormInt") || ! dspec["ParticleANormInt"].isGroup())
	throw RLException("Could not look up ParticleANormInt.");
  Setting &partA = dspec["ParticleANormInt"];
  if( ! dspec.exists("ParticleBNormInt") || ! dspec["ParticleBNormInt"].isGroup())
	throw RLException("Could not look up ParticleBNormInt.");
  Setting &partB = dspec["ParticleBNormInt"];

  if( ! dspec.lookupValue("LegendreRule", useLegendreRule) )
	throw RLException("Could not look up UseLegendreRule.");
  
  
  if(! partA.lookupValue("Start", particleDomain.at(0).start) )
	throw RLException("Could not look up time start.");
  if(! partA.lookupValue("Stop", particleDomain.at(0).stop))
	throw RLException("Could not look up time stop.");
  SizetLookup(partA, "Precision", particleDomain.at(0).precision);
  
  if(! partB.lookupValue("Start", particleDomain.at(1).start) )
	throw RLException("Could not look up time start.");
  if(! partB.lookupValue("Stop", particleDomain.at(1).stop))
	throw RLException("Could not look up time stop.");
  SizetLookup(partB, "Precision", particleDomain.at(1).precision);
}

void FluxConfig::InitUnits(Setting & root)
{
  Setting & units = root["Units"];
  if( ! units.lookupValue("HbarTimesLambda", myUnits.hbarTimesLambda) )
	{
	  throw RLException("HbarTimesLambda unit not properly specified.");
	}
  if( ! units.lookupValue("MassOverLambda2", myUnits.massOverLambda2) )
	{
	  throw RLException("MassOverLambda2 unit not properly specified.");
	}  
}

FluxConfig::~FluxConfig()
{
  
}


uint FluxConfig::GetVerbosityLevel() const
{
  return verbosityLevel;
}

uint FluxConfig::GetNumberOfThreads() const
{
  return numberOfThreads;
}

const InputFiles & FluxConfig::GetInputFiles() const
{
  return myInputFiles;
}

const vector<DomainSpecific> & FluxConfig::GetParticleDomain() const
{
  return particleDomain;
}

const UnitSet & FluxConfig::GetUnits() const
{
  return myUnits;
}

const vector<size_t> & FluxConfig::GetSingleParticleBasisNumber() const
{
  return mySingleParticleBasisNumber;
}

const string & FluxConfig::GetOutputDensityFile() const
{
  return outputDensity;
}

const string & FluxConfig::GetOutputFluxFile() const
{
  return outputFlux;
}

const bool FluxConfig::GetUseLegendreRule() const
{
  return useLegendreRule;
}
