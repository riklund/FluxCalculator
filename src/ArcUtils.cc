#include "ArcUtils.hh"

double InterpolateFlux(const vector<double> & xvals, const vector<vector<double> > fVals, double x, double y)
{
  double delta = xvals.at(1) - xvals.at(0);
  size_t xUp = distance(xvals.begin(), upper_bound(xvals.begin(), xvals.end(), x));
  size_t yUp = distance(xvals.begin(), upper_bound(xvals.begin(), xvals.end(), y));
  if(xUp == 0 || yUp == 0)
	throw RLException("xLow or yLow was zero.");
  size_t xLo = xUp - 1;
  size_t yLo = yUp - 1;


  if(xUp == xvals.size())
	--xUp;
  if(yUp == xvals.size())
	--yUp;

  ///Interpolate as a plane. 

  double dx = x - xvals.at(xLo);
  double dy = y - xvals.at(yLo);

  double fx0y0 = fVals.at(xLo).at(yLo);
  double dFdx = (fVals.at(xUp).at(yLo) - fVals.at(xLo).at(yLo)) / delta;
  double dFdy = (fVals.at(xLo).at(yUp) - fVals.at(xLo).at(yLo)) / delta;


  double f = fx0y0 + dFdx * dx + dFdy * dy;
  return f;
}


double TrapezoidIntegrate(const vector<double> & vect, double delta)
{
  double sum = 0.0;
  for(size_t i = 0; i<vect.size(); ++i)
	{
	  double preFact = delta;
	  if(i == 0 || i == vect.size() - 1)
		preFact *= 0.5;
	  sum += preFact * vect.at(i);
	}
  return sum;
}

double TrapezoidIntegrate(const vector<vector<double> > & vect, double delta )
{
  vector<double> interm(vect.size());
  for(size_t i = 0; i<vect.size(); ++i)
	{
	  interm.at(i) = TrapezoidIntegrate(vect.at(i), delta);
	}

  return TrapezoidIntegrate(interm, delta);
}


double ArcIntegrate(const vector<vector<double> > & vect, const vector<double> & xvals)
{
  if(vect.size() != xvals.size())
	throw RLException("vect-size (%d) xvals-size (%d) mismatch.", vect.size(), xvals.size());

  double delta = xvals.at(1) - xvals.at(0);
  double R = xvals.back() - xvals.front();
  vector<double> interm(vect.size());
  for(size_t i = 0; i<vect.size(); ++i)
	{
	  interm.at(i) = 0.0;
	  double x = xvals.at(i); 
	  double ymax = sqrt(pow(R, 2) - pow(x,2));
	  size_t idy = distance(xvals.begin(), upper_bound(xvals.begin(), xvals.end(), ymax-.5*delta));
	  for(size_t j = 0; j<idy; ++j)
		{
		  double preFact = delta;
		  if(j==0 || j == idy - 1)
			preFact *= 0.5;
		  interm.at(i) += preFact * vect.at(i).at(j);
		}
	}
  return TrapezoidIntegrate(interm, delta);
}

