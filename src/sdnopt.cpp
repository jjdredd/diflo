/*****************************************************************************
 * This program is free software: you can redistribute it and/or modify	     *
 * it under the terms of the GNU General Public License version 2.	     *
 * 									     *
 * Author: Rahim Usubov <usubov@theor.jinr.ru>				     *
 * 									     *
 * GNU GSL is required to build this program. Build details are provided in  *
 * the Makefile.							     *
 *****************************************************************************/

#include "handedness.hpp"
#include "data.hpp"

#include <cmath>
#include <cunistd>
#include <gsl/gsl_statistics_double.h>

#include <fstream>
#include <string>
#include <sstream>
/* HSD and PHSD seem to have slightly different output formats,
 * to take this into account, a PHSD preprocessor variable is used */
/* #define PHSD */

void print_usage() {
	printf ("");
	return;
}

int main(int argc, char** argv){
  
	char c, *file_data = NULL, *file_input = NULL;
	DataVersion dversion;
	bool alice = false, req_input = true;
	while ((c = getopt (argc, argv, "A:P:H:i:")) != -1) {

		switch (c) {

		case 'A':	// ALICE data
			file_data = optarg;
			alice = true;
			req_input = false;
			break;

		case 'P':	// PHSD data
			dversion = HSD_VER_PHSD;
			file_data = optarg;
			break;

		case 'H':	// HSD original data
			dversion = HSD_VER_ORIG;
			file_data = optarg;
			break;

		case 'i':
			file_input = optarg;

		case '?':
		default:
			print_usage();
			return -1;
		}
	}

	if (req_input && !alice) {
		print_usage();
		return -1;
	}

	if (alice) {
		std::ofstream files[8];
		for (unsigned i = 0; i < 8; i++) {
			std::string name = "res/oct_multip_"
				+ std::to_string(i);
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
				mean = gsl_stats_mean (&etas[i][0], 1,
						       numevents);
				rsdm = rsn * gsl_stats_sd_m (&etas[i][0], 1,
							     numevents, mean);
				files[i] << rsn << '\t' << start
					 << '\t' << mean
					 << '\t' << rsdm << std::endl;
			}
		}

	} else {
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
		std::vector<double> etas[2], evet;
		std::vector<unsigned> evnm;

		etas[0].reserve(D.ISUBS * D.NUM);
		etas[1].reserve(D.ISUBS * D.NUM);
		evet.reserve(2);
		evnm.reserve(2);

#if 0
		for(unsigned isub = 0; isub < D.ISUBS; isub++){
			for(unsigned irun = 0; irun < D.NUM; irun++){
				EventEta(D.P[isub][irun], evet, evnm, rpa);
				etas[0].push_back(evet[0]);
				etas[1].push_back(evet[1]);
			}
		}

		for(unsigned i = 0; i < 2; i++){
			double rsn, mean, rsdm;
			unsigned numevents = D.ISUBS * D.NUM;
			rsn = 1/sqrt(numevents);
			mean = gsl_stats_mean (&etas[i][0], 1, numevents);
			rsdm = rsn * gsl_stats_sd_m (&etas[i][0], 1,
						     numevents, mean);
			std::cout << rsn << '\t' << mean
				  << '\t' << rsdm << std::endl;
		}
		std::cout << "Correlation: "
			  << gsl_stats_correlation (&etas[0][0], 1,
						    &etas[1][0], 1,
						    D.NUM * D.ISUBS)
			  << std::endl;
#endif	// end of old calc code

		std::ofstream fout("res_vs_rpa.txt");
		for(double rpa = 0; rpa < M_PI; rpa += 0.4) {
			etas[0].clear();
			etas[1].clear();
			etas[0].reserve(D.ISUBS * D.NUM);
			etas[1].reserve(D.ISUBS * D.NUM);
			for(unsigned isub = 0; isub < D.ISUBS; isub++){
				for(unsigned irun = 0; irun < D.NUM; irun++){
					EventEta(D.P[isub][irun], evet,
						 evnm, rpa);
					etas[0].push_back(evet[0]);
					etas[1].push_back(evet[1]);
				}
			}

			fout << rpa << '\t';
			for(unsigned i = 0; i < 2; i++){
				double rsn, mean, rsdm;
				unsigned numevents = D.ISUBS * D.NUM;
				rsn = 1/sqrt(numevents);
				mean = gsl_stats_mean (&etas[i][0], 1,
						       numevents);
				rsdm = rsn * gsl_stats_sd_m (&etas[i][0], 1,
							     numevents, mean);
				fout << mean << '\t' << rsdm << '\t';
			}
			fout << gsl_stats_correlation (&etas[0][0], 1,
						       &etas[1][0], 1,
						       D.NUM * D.ISUBS)
			     << std::endl;
		}

	}

	return 0;
}
