#include "TestClassTest.hh"

int TestClassTest::TestCase1() const
{
  return 0; ///Normally, perform tests here. Return nonzero if required.
}

int TestClassTest::RunTests() const
{
  int code1 = TestCase1();
  if(code1)
    return code1; ///Return nonzero for incorrect behavior (test failure).

  return 0; ///Return 0 for correct behavior.
}


string TestClassTest::ToString() const
{
  return "Dummy test.";
}
