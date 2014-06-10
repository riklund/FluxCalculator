#include "ArcUtilsTest.hh"

bool ArcUtilsTest::InterpolateFlux_PerformsCorrectly() const
{
  vector<double> xvals;
  vector<vector<double> > fvals(3,vector<double>(3, 0.0));
  xvals.push_back(1);
  xvals.push_back(2);
  xvals.push_back(3);
  for(int i = 0; i<9; ++i)
	{
	  fvals[i/3][i%3] = i;
	}

  for(int i = 0; i<9; ++i)
	{
	  if(!DBL_EQUAL(InterpolateFlux(xvals, fvals, 1+i/3, 1+i%3), i) )
		return false;
	}

  if(!DBL_EQUAL(InterpolateFlux(xvals, fvals, 1.5, 1.5), 2.0))
	return false;
  return true;
}


bool ArcUtilsTest::TrapezoidIntegrate_PerformsCorrectly() const
{
  vector<double> yvals;
  double delta = 1;
  for(size_t i = 0; i<11; ++i)
	yvals.push_back(1);

  if(TrapezoidIntegrate(yvals, delta) != 10)
	return false;

  vector<vector<double> > dyVals(11, yvals);

  if(TrapezoidIntegrate(dyVals, delta) != 100)
	return false;

  return true;
}

bool ArcUtilsTest::ArcIntegrate_PerformsCorrectly() const
{
  size_t ulim = 10001;

  vector<double> xvals(ulim);
  vector<vector<double> > fvals(ulim, vector<double>(ulim,0.0));
  for(size_t i = 0; i<ulim; ++i)
	{
	  xvals.at(i) = ((double)i)/(ulim -1);
	  for(size_t j = 0; j<ulim; ++j)
		{
		  fvals.at(i).at(j) = 1.0;
		}
	}
  if(abs(ArcIntegrate(fvals, xvals) - PI/4) > 1E-2)
	return false;

  return true;
}


int ArcUtilsTest::RunTests() const
{
  if(!InterpolateFlux_PerformsCorrectly())
	return 1;
  if(!TrapezoidIntegrate_PerformsCorrectly()) 
	return 2;
  if(!ArcIntegrate_PerformsCorrectly())
	return 3;


  return 0; ///Return 0 for correct behavior.
}


string ArcUtilsTest::ToString() const
{
  return "ArcUtils";
}
