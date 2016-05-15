#include "handedness.hpp"
#include "dist.hpp"
#include "data.hpp"

#include <cmath>
#include <unistd.h>
#include <gsl/gsl_statistics_double.h>

#include <fstream>
#include <string>
#include <sstream>

int main(int argc, char** argv){
  

	if (argc < 4) {
		std::cerr << "too few arguments" << std::endl;
		return -1;
	}

	std::ifstream s(argv[1]);
	DataHSD D(s);
	s.close();

	// mesons
	{
		s.open(argv[2]);
		if(s.is_open()) D.readin_particles(s, true);
		else std::cout << "Warning: couldn't open mesons file "
			       << argv[2] << '\n';
		s.close();
	}

	// baryons
	{
		s.open(argv[3]);
		if(s.is_open()) D.readin_particles(s, false);
		else std::cout << "Warning: couldn't open baryons file "
			       << argv[3] << '\n';
		s.close();
	}

	std::vector<double> rpa;
	for(unsigned isub = 0; isub < D.ISUBS; isub++){
		for(unsigned irun = 0; irun < D.NUM; irun++){
			// double angle = RPA_by_multip(D.P[isub][irun]);
			double angle = RPA_by_multip(D.P[isub][irun]);
			rpa.push_back(angle);
			std::cout << angle << std::endl;
		}
	}

	double rsn, mean, rsdm;
	unsigned numevents = D.ISUBS * D.NUM;
	rsn = 1/sqrt(numevents);
	mean = gsl_stats_mean (&rpa[0], 1, numevents);
	rsdm = rsn * gsl_stats_sd_m (&rpa[0], 1, numevents, mean);
	std::cout << rsn << '\t' << mean
		  << '\t' << rsdm << std::endl;
	return 0;
}
