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
#include "dist.hpp"
#include "data.hpp"

#include <cmath>
#include <unistd.h>
#include <gsl/gsl_statistics_double.h>

#include <fstream>
#include <string>
#include <sstream>

void print_usage() {
	printf ("sdnopt.elf: calculate handedness for various data\n"
		"\t -A <file> - parse ALICE data from <file>\n\n"
		"\t -P <file> - parse PHSD data from <file>\n\n"
		"\t -H <file> - parse HSD data from <file>\n\n"
		"The following options are only for (P)HSD data.\n\n"
		"\t -i <file> - use <file> as input file for"
		"(P)HSD program\n\n"
		"\t -t <int> - use only particle type <int>\n\n"
		"\t -m - ALICE angle analysis\n\n"
		"\t -2 - calculate handedness in diants\n\n"
		"\t -8 - calculate handedness for octants\n\n"
		"\t -a - calculate handedness vs angle (implies -2)\n\n");
	return;
}

int main(int argc, char** argv){
  
	char c, *file_data = NULL, *file_input = NULL;
	bool alice = false, req_input = true, pick = false,
		angle = false, in_octants = true, angle_anal = false;
	DataVersion dversion;
	int type = 0;
	while ((c = getopt (argc, argv, "A:P:H:i:t:m28a")) != -1) {

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
			break;

		case 't':
			pick = true;
			type = std::stoi(optarg);
			break;

		case '2':
			in_octants = false;
			angle = false;
			break;

		case '8':
			in_octants = true;
			angle = false;
			break;

		case 'a':
			angle = true;
			in_octants = false;
			break;

		case 'm':
			angle_anal = true;
			break;

		case '?':
		default:
			print_usage();
			return -1;
		}
	}

	if (req_input && !alice && !file_input) {
		print_usage();
		return -1;
	}

	if (alice && angle_anal) {

		unsigned w = 50;
		std::string base_fname = "res/aaa.txt_";
		std::ofstream dep_out;
		dep_out.open("res/ratio_vs_n.txt", std::ofstream::out);

		for (unsigned start = 250; start < 650; start += w) {
			std::ofstream ofile;
			ofile.open(base_fname + std::to_string(start),
				   std::ofstream::out);

			std::vector<double> RatioS, RPAngleS;
			ALICEData D(file_data);
			event e;

			std::cout << "start " << start << ", width "
				  << w << ", output file "
				  << base_fname + std::to_string(start)
				  << std::endl;
			while (D.FetchNumEvent(e, start, start + w)) {
				double ratio, max_angle, rpangle;
				ratio = MaxHandedRatio(e, max_angle);
				rpangle = RPA_by_multip(e);
				RatioS.push_back(ratio);
				RPAngleS.push_back(max_angle - rpangle);
				std::cout << '\t' << "ratio " << ratio
					  << ", diff "
					  << max_angle - rpangle
					  << std::endl;
			}
			for (unsigned i = 0; i < RatioS.size(); i++)
				ofile << RatioS[i] << '\t'
				      << RPAngleS[i] << std::endl;

			double mean, rsn, rsdm;
			rsn = 1 / sqrt(start);
			mean = gsl_stats_mean (&RatioS[0], 1, RatioS.size());
			rsdm = rsn * gsl_stats_sd_m (&RatioS[0], 1,
						     RatioS.size(), mean);
			dep_out << rsn << '\t' << start << '\t'
				<< mean << '\t' << rsdm << std::endl;

			ofile.close();
		}

		return 0;
	}


	// ALICE doesn't work!!!
	// needs a fix!!
	// XXX: FIXME
	if (alice) {

		HandednessExp H;
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

			ALICEData D(file_data);
			event e;

			while (D.FetchNumEvent(e, start, start + w)) {
				H.RPAngle = 0;
				H.EventEta(e, evet, evnm);
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

		return 0;
	}

	std::ifstream s(file_input);
	DataHSD *D = NULL;

	//
	// DON'T FORGET ABOUT PARTICLE TYPE HERE!
	//

	switch (dversion) {

		// D is leaking everywhere!

	case HSD_VER_ORIG:
		D = new DataHSD(s, pick, type);
		break;

	case HSD_VER_COORD:
		D = new DataHSDC(s, pick, type);
		break;

	case HSD_VER_PHSD:
		D = new DataPHSD(s, pick, type);
		break;

	}
	s.close();

	// mesons
	{
		s.open(file_data);
		if(s.is_open()) D->readin_particles(s, true);
		else std::cout << "Warning: couldn't open mesons file "
			       << file_data << '\n';
		s.close();
	}

	if (in_octants) {

		std::vector<double> etas[8], evet;
		std::vector<unsigned> evnm;
		Handedness H;

		for (unsigned oct = 0; oct < 8; oct++)
			etas[oct].reserve(D->ISUBS * D->NUM);

		evet.reserve(8);
		evnm.reserve(8);

		for(unsigned isub = 0; isub < D->ISUBS; isub++){
			for(unsigned irun = 0; irun < D->NUM; irun++){
				H.EventEta(D->P[isub][irun], evet, evnm);

				for (unsigned oct = 0; oct < 8; oct++)
					etas[oct].push_back(evet[oct]);
			}
		}

		for(unsigned i = 0; i < 8; i++){
			double rsn, mean, rsdm;
			unsigned numevents = D->ISUBS * D->NUM;
			rsn = 1/sqrt(numevents);
			mean = gsl_stats_mean (&etas[i][0], 1, numevents);
			rsdm = rsn * gsl_stats_sd_m (&etas[i][0], 1,
						     numevents, mean);
			std::cout << rsn << '\t' << mean
				  << '\t' << rsdm << std::endl;
		}

		// forget about correlations for a while
		// std::cout << "Correlation: "
		// 	  << gsl_stats_correlation (&etas[0][0], 1,
		// 				    &etas[1][0], 1,
		// 				    D->NUM * D->ISUBS)
		// 	  << std::endl;

		return 0;
	}

	std::vector<double> etas[2], evet;
	std::vector<unsigned> evnm;
	HandednessExp H;

	evet.reserve(2);
	evnm.reserve(2);

	if (angle) {

		std::ofstream fout("res_vs_rpa.txt");
		for(double rpa = 0; rpa < M_PI; rpa += 0.4) {
			etas[0].clear();
			etas[1].clear();
			etas[0].reserve(D->ISUBS * D->NUM);
			etas[1].reserve(D->ISUBS * D->NUM);
			for(unsigned isub = 0; isub < D->ISUBS; isub++){
				for(unsigned irun = 0; irun < D->NUM; irun++){
					H.RPAngle = rpa;
					H.EventEta(D->P[isub][irun], evet,
						    evnm);
					etas[0].push_back(evet[0]);
					etas[1].push_back(evet[1]);
				}
			}

			fout << rpa << '\t';
			for(unsigned i = 0; i < 2; i++){
				double rsn, mean, rsdm;
				unsigned numevents = D->ISUBS * D->NUM;
				rsn = 1/sqrt(numevents);
				mean = gsl_stats_mean (&etas[i][0], 1,
						       numevents);
				rsdm = rsn * gsl_stats_sd_m (&etas[i][0], 1,
							     numevents, mean);
				fout << mean << '\t' << rsdm << '\t';
			}
			fout << gsl_stats_correlation (&etas[0][0], 1,
						       &etas[1][0], 1,
						       D->NUM * D->ISUBS)
			     << std::endl;
		}

	} else {

		for(unsigned isub = 0; isub < D->ISUBS; isub++){
			for(unsigned irun = 0; irun < D->NUM; irun++){
				H.RPAngle = 0;
				H.EventEta(D->P[isub][irun], evet, evnm);
				etas[0].push_back(evet[0]);
				etas[1].push_back(evet[1]);
			}
		}

		for(unsigned i = 0; i < 2; i++){
			double rsn, mean, rsdm;
			unsigned numevents = D->ISUBS * D->NUM;
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
						    D->NUM * D->ISUBS)
			  << std::endl;
	}

	delete D;
	return 0;
}
