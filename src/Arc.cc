#include "Arc.hh"

int main(int argc, char ** argv)
{
  if(argc != 7)
	{
	  cerr << "Usage: " << argv[0] << " input_flux input_density a b N output_file" << endl;
	  return 1;
	}
  double a = atof(argv[3]);
  double b = atof(argv[4]);
  size_t N = atol(argv[5]);
  double delta = (b-a)/N;

  VerbosePrinter myPrinter(100);
  vector<double> xvals(N+1);
  for(size_t i = 0; i<N+1; ++i)
	{
	  xvals.at(i) = a + delta * i;
	}


  myPrinter.Print(2, "Obtaining density ratio...");
  double densRatio = GetDensityRatio(xvals, argv[2]);
  myPrinter.Print(2, "done.\n");
  myPrinter.Print(1, "Density ratio is %lf.\n", densRatio);

  myPrinter.Print(2, "Obtaining and writing arc...");
  WriteArcToFile(xvals, argv[1], argv[6], densRatio);
  myPrinter.Print(2, "done.\n");

  return 0;
}

void WriteArcToFile(const vector<double> & xvals, const char * inputFileName, const char * outputFileName, double densRatio)
{
  vector<vector<double> > fVals(xvals.size(), vector<double>(xvals.size(), INIT_VAL));
  FILE * fin = fopen(inputFileName, "r");
  if(fin == NULL)
	throw RLException("Could not open file '%s'", inputFileName);
  double x1, x2, fluxX1, fluxX2;

  while(fscanf(fin, "%lf %lf %lf %lf ", &x1, &x2, &fluxX1, &fluxX2) != EOF)
	{
	  size_t idx1 = distance(xvals.begin(), upper_bound(xvals.begin(), xvals.end(), x1-EPS));
	  size_t idx2 = distance(xvals.begin(), upper_bound(xvals.begin(), xvals.end(), x2-EPS));
	  if(idx1 >= xvals.size() || idx2 >= xvals.size())
		{
		  throw RLException("Invalid read : %13.10e %13.10e %13.10e %13.10e", x1, x2, fluxX1, fluxX2);
		}
	  if(abs(fVals.at(idx1).at(idx2) - INIT_VAL) > EPS)
		{
		  throw RLException("Visited same value twice: %lf, %lf, was %+13.10e, attempted %+13.10e.", x1, x2, fVals.at(idx1).at(idx2), fluxX1 + fluxX2);
		}
	  fVals.at(idx1).at(idx2) = fluxX1 + fluxX2;
	}
  fclose(fin);

  double dTheta = PI * 0.5 / (NUMBER_OF_THETA);
  double R = xvals.back() - xvals.front();
  FILE * fout = fopen(outputFileName, "w");
  if(fout == NULL)
	throw RLException("Could not open file '%s'", outputFileName);

  for(size_t i = 0; i<NUMBER_OF_THETA; ++i)
	{
	  double theta = dTheta * i;
	  double x = xvals.at(0) + R * cos(theta);
	  double y = xvals.at(0) + R * sin(theta);

	  double flux = InterpolateFlux(xvals, fVals, x, y);
	  fprintf(fout, "%+13.10e %+13.10e\n", theta, flux * densRatio);
	}
  fclose(fout);
}


double GetDensityRatio(const vector<double> & xvals, const char * fileName)
{
  vector<vector<double> > dens(xvals.size(), vector<double>(xvals.size(), INIT_VAL));

  ///Read density from file.
  double x1, x2, dval;
  FILE * densityFile = fopen(fileName, "r");
  if(densityFile == NULL)
	throw RLException("Could not open file '%s'", densityFile);

  while(fscanf(densityFile, "%lf %lf %lf ", &x1, &x2, &dval) != EOF)
	{
	  size_t idx1 = distance(xvals.begin(), upper_bound(xvals.begin(), xvals.end(), x1-EPS));
	  size_t idx2 = distance(xvals.begin(), upper_bound(xvals.begin(), xvals.end(), x2-EPS));
	  if(idx1 >= xvals.size() || idx2 >= xvals.size())
		{
		  throw RLException("Invalid read : %13.10e %13.10e %13.10e", x1, x2, dval);
		}
	  if(abs(dens.at(idx1).at(idx2) - INIT_VAL) > EPS)
		throw RLException("Visited same value twice: %lf, %lf, was %+13.10e, attempted %+13.10e.", x1, x2, dens.at(idx1).at(idx2), dval);
	  dens.at(idx1).at(idx2) = dval;
	}
  fclose(densityFile);

  ///Integrate. We are using the trapetzoid rule here.
  double fullArea = TrapezoidIntegrate(dens, xvals.at(1) - xvals.at(0));
  double arcArea = ArcIntegrate(dens, xvals);
  return fullArea / arcArea;
}
