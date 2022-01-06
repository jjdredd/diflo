#include "statistics.hpp"

#include <cmath>
#include <unistd.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

const unsigned SZAR = 100;

int main(int argc, char** argv){
  
	double *a, *b;
	double gslCor, customCor;
	StatAggregatorCorrelations<double> csac(0);
	gsl_rng *rng = gsl_rng_alloc(gsl_rng_taus);
	gsl_rng_set(rng, 1672);

	a = new double[SZAR];
	b = new double[SZAR];

	// test goes here

	// fill a, b
	for(unsigned i = 0; i < SZAR; i++) {
		b[i] = gsl_ran_gaussian(rng, 0.5);
		a[i] = gsl_ran_exponential(rng, 0.25);
	}

	// gsl correlation
	gslCor = gsl_stats_correlation(a, 1, b, 1, SZAR);

	// custom correlation
	for(unsigned i = 0; i < SZAR; i++) {
		csac.ConsumeDataPoint(a[i], b[i]);
	}
	customCor = csac.CurrentValue();

	std::cout << "gsl correlation is: \t" << gslCor << std::endl
		  << "statistics.hpp correlations is: \t" << customCor << std::endl
		  << std::endl;

	delete[] a;
	delete[] b;
	gsl_rng_free(rng);

	return 0;
}
