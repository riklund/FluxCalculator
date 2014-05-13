#include "ManagedOutputFiles.hh"

ManagedOutputFiles::ManagedOutputFiles(VerbosePrinter & myPrinter, const OutputFiles & data)
  :TimeSeries(NULL), SingleBasis(2, NULL)
{
  for(uint i = 0; i<2; ++i)
	{
	  VerboseOpen(myPrinter, &SingleBasis.at(i), data.SingleParticle.at(i), "SingleParticle");
	}
  VerboseOpen(myPrinter, &TimeSeries, data.TimeSeries, "TimeSeries");
}

ManagedOutputFiles::~ManagedOutputFiles()
{
  if(TimeSeries)
	fclose(TimeSeries);
  for(vector<FILE*>::const_iterator it = SingleBasis.begin(); it!=SingleBasis.end(); ++it)
	{
	  if(*it)
		fclose(*it);
	}
}

void ManagedOutputFiles::VerboseOpen(VerbosePrinter & myPrinter, FILE ** fout, string fileName, const char * description)
{
  if(fileName.empty())
	{
	  myPrinter.Print(4, "Empty filename, will not save '%s'.\n", description);
	}
  else
	{
	  myPrinter.Print(4, "Will save '%s' to file '%s'.\n", description, fileName.c_str());
	  (*fout) = fopen(fileName.c_str(), "w");
	  if((*fout) == NULL)
		{
		  throw RLException("Could not open output file '%s' for '%s'.", fileName.c_str(), description);
		}
	}
}
