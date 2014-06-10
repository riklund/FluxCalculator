#ifndef ArcUtilsTest_hh
#define ArcUtilsTest_hh 1

///This should always be included
#include "GenericUnitTest.hh"

///General inclusion statements.
#include <iostream>
#include <assert.h>
#include <map>
#include <vector>
#include <cmath>
#include <string.h>
#include <stdio.h>
#include <algorithm>

#include "ArcUtils.hh"

#ifndef DBL_EQUAL
#define DBL_EQUAL(a, b) (abs(b-a) < 1E-9)
#endif

using namespace std;

/**
   Test class.
 */
//It is imperative to have this on a single line before the inheritance line.
class ArcUtilsTest 
  : public GenericUnitTest ///Always extend GenericUnitTest.
{
 public:
  int RunTests() const; ///Main function.
  string ToString() const; ///Should return a descriptive name of the function.
 protected:
  bool InterpolateFlux_PerformsCorrectly() const;
  bool TrapezoidIntegrate_PerformsCorrectly() const;
  bool ArcIntegrate_PerformsCorrectly() const;
};
#endif
