#ifndef ManagedOutputFiles_hh
#define ManagedOutputFiles_hh 1

#include "VerbosePrinter.hh"
#include "DensityConfig.hh"
#include <cstdio>
#include <vector>

using namespace std;

class ManagedOutputFiles
{
public:
  ManagedOutputFiles(VerbosePrinter & myPrinter, const OutputFiles & data);
  ~ManagedOutputFiles();
protected:
  static void VerboseOpen(VerbosePrinter & myPrinter, FILE ** fout, string fileName, const char * description);
  
public:
  FILE * TimeSeries;
  vector<FILE*> SingleBasis;
};

#endif
