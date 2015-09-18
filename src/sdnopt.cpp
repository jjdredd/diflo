/*****************************************************************************
 * This program is free software: you can redistribute it and/or modify	     *
 * it under the terms of the GNU General Public License version 2.	     *
 * 									     *
 * Author: Rahim Usubov <usubov@theor.jinr.ru>				     *
 * 									     *
 * GNU GSL is required to build this program. Build detailes are provided in *
 * the Makefile.							     *
 *****************************************************************************/

#include "handedness.hpp"
#include "data.hpp"

#include <cmath>
#include <gsl/gsl_statistics_double.h>

#include <fstream>
#include <string>
#include <sstream>
/* HSD and PHSD seem to have slightly different output formats,
 * to take this into account, a PHSD preprocessor variable is used */
/* #define PHSD */

int main(int argc, char** argv){
  
	if (argc < 2){
		printf ("not enough arguments\n"
			"usage:\n"
			"first arg input file\n");
		return -1;
	}

	std::ofstream files[8];
	for (unsigned i = 0; i < 8; i++) {
		std::string name = "res/oct_multip_" + std::to_string(i);
		files[i].open(name, std::ofstream::out);
	}

	unsigned w = 50;
	for (unsigned start = 50; start < 600; start += w) {
		std::vector<double> etas[8], evet;
		std::vector<unsigned> evnm;

		ALICEData D(argv[1]);
		event e;

		while (D.FetchNumEvent(e, start, start + w)) {
			EventEta(e, evet, evnm);
			for(unsigned i = 0; i < 8; i++) {
				etas[i].push_back(evet[i]);
			}
		}

		for(unsigned i = 0; i < 8; i++){
			double rsn, mean, rsdm;
			unsigned numevents = etas[i].size();
			rsn = 1/sqrt(numevents);
			mean = gsl_stats_mean (&etas[i][0], 1, numevents);
			rsdm = rsn * gsl_stats_sd_m (&etas[i][0], 1,
						     numevents, mean);
			files[i] << rsn << '\t' << start
				 << '\t' << mean
				 << '\t' << rsdm << std::endl;
		}
	}

#if 0				// this is for HSD
	std::ifstream s(argv[2]);
	data D(s, HSD_VER_ORIG, true, 1);
	s.close();
	// mesons
	{
		s.open(argv[1]);
		if(s.is_open()) D.readin_particles(s, true);
		// else std::cout << "Warning: couldn't open mesons file "
		// 		   << argv[1] << '\n';
		s.close();
	}
	for(unsigned i = 0; i < 8; i++)
		etas[i].reserve(D.ISUBS * D.NUM);
	evet.reserve(8);
	evnm.reserve(8);
	for(unsigned isub = 0; isub < D.ISUBS; isub++){
		for(unsigned irun = 0; irun < D.NUM; irun++){
			EventEta(D.P[isub][irun], evet, evnm);
			for(unsigned i = 0; i < 8; i++)
				etas[i].push_back(evet[i]);
		}
	}

	for(unsigned i = 0; i < 8; i++){
		double rsn, mean, rsdm;
		unsigned numevents = D.ISUBS * D.NUM;
		rsn = 1/sqrt(numevents);
		mean = gsl_stats_mean (&etas[i][0], 1, numevents);
		rsdm = rsn * gsl_stats_sd_m (&etas[i][0], 1,
					     numevents, mean);
		std::cout << rsn << '\t' << 1.0/numevents
			  << '\t' << mean << '\t' << fabs(mean)
			  << '\t' << rsdm << std::endl;
	}
#endif	// end of testing code for HSD data

	return 0;
}
