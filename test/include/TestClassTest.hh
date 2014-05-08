#ifndef TestClassTest_hh
#define TestClassTest_hh 1

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

using namespace std;

/**
   Test class.
 */
//It is imperative to have this on a single line before the inheritance line.
class TestClassTest 
  : public GenericUnitTest ///Always extend GenericUnitTest.
{
 public:
  int RunTests() const; ///Main function.
  string ToString() const; ///Should return a descriptive name of the function.
 protected:
  int TestCase1() const; /// Test case.
};
#endif
