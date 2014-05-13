#include "DensityConfig.hh"

DensityConfig::DensityConfig()
  :verbosityLevel(0), numberOfThreads(1), particleDomain(2)
{
  
}

DensityConfig::DensityConfig(const char * fileName)
 : particleDomain(2)
{
  ReadFile(fileName);
}

void DensityConfig::ReadFile(const char * fileName)
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
  InitOutputFiles(root);
  InitInputFiles(root);
  InitDomainSpecifics(root);
  InitUnits(root);
}

void DensityConfig::InitProgramGenerals(Setting & root)
{
  if( ! root.exists("Program") || !root["Program"].isGroup())
	throw RLException("Could not look up group 'Program' in config file.");
  Setting & program = root["Program"];
  if( ! program.lookupValue("VerbosityLevel", verbosityLevel) )
	throw RLException("Could not look up verbosity level.");
  if( ! program.lookupValue("NumberOfThreads", numberOfThreads) )
	throw RLException("Could not look up number of threads.");
}

void DensityConfig::InitInputFiles(Setting & root)
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

void DensityConfig::InitOutputFiles(Setting & root)
{
  if( ! root.exists("OutputFiles") || ! root["OutputFiles"].isGroup() )
	throw RLException("Invalid output file group");
  Setting & output = root["OutputFiles"];
  if( ! output.lookupValue("TimeSeries", myOutputFiles.TimeSeries) )
	throw RLException("Could not look up value for TimeSeries output file.");
  if( ! output.lookupValue("FullFunctions", myOutputFiles.FullFunctions) )
	throw RLException("Could not look up value for FullFunctions output file.");
  if( ! output.lookupValue("SingleA", myOutputFiles.SingleParticle.at(0)) )
	throw RLException("Could not look up value for SingleA output file.");
  if( ! output.lookupValue("SingleB", myOutputFiles.SingleParticle.at(1)) )
	throw RLException("Could not look up value for SingleB output file.");
  if( ! output.lookupValue("SingleABasis", myOutputFiles.SingleBasis.at(0)))
	throw RLException("Could not look up value for SingleABasis output property.");
  if( ! output.lookupValue("SingleBBasis", myOutputFiles.SingleBasis.at(1)))
	throw RLException("Could not look up value for SingleBBasis output property.");
}

void DensityConfig::SizetLookup(const Setting & root, const char * propName, size_t & output)
{
  unsigned int tmp;
  if(!root.lookupValue(propName, tmp))
	{
	  throw RLException("Lookup of property '%s' failed.", propName);
	}
  output=tmp;
  return;
}

void DensityConfig::InitDomainSpecifics(Setting & root)
{
    if( ! root.exists("DomainSpecifics") || ! root["DomainSpecifics"].isGroup() )
	  throw RLException("Invalid DomainSpecifics file group");
	Setting & dspec = root["DomainSpecifics"];
	if( ! dspec.exists("Time") || ! dspec["Time"].isGroup())
	  throw RLException("Could not look up Time.");
	Setting & time = dspec["Time"];
	if( ! dspec.exists("ParticleA") || ! dspec["ParticleA"].isGroup())
	  throw RLException("Could not look up ParticleA.");
	Setting &partA = dspec["ParticleA"];
	if( ! dspec.exists("ParticleB") || ! dspec["ParticleB"].isGroup())
	  throw RLException("Could not look up ParticleB.");
	Setting &partB = dspec["ParticleB"];


	if(! time.lookupValue("Start", timeDomain.start) )
	  throw RLException("Could not look up time start.");
	if(! time.lookupValue("Stop", timeDomain.stop))
	  throw RLException("Could not look up time stop.");
	SizetLookup(time, "Precision", timeDomain.precision);
	
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

void DensityConfig::InitUnits(Setting & root)
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

DensityConfig::~DensityConfig()
{
  
}


uint DensityConfig::GetVerbosityLevel() const
{
  return verbosityLevel;
}

uint DensityConfig::GetNumberOfThreads() const
{
  return numberOfThreads;
}

const InputFiles & DensityConfig::GetInputFiles() const
{
  return myInputFiles;
}

const OutputFiles & DensityConfig::GetOutputFiles() const
{
  return myOutputFiles;
}

const vector<DomainSpecific> & DensityConfig::GetParticleDomain() const
{
  return particleDomain;
}

const DomainSpecific & DensityConfig::GetTimeDomain() const
{
  return timeDomain;
}

const UnitSet & DensityConfig::GetUnits() const
{
  return myUnits;
}

