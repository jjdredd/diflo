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
/* HSD and PHSD seem to have slightly different output formats,
 * to take this into account, a PHSD preprocessor variable is used */
/* #define PHSD */

int main(int argc, char** argv){
  
	if (argc < 3){
		printf ("not enough arguments\n"
			"usage:\n"
			"first arg fort.301\n"
			"second arg input\n");
		return -1;
	}

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

	std::vector<double> etas[8], evet;
	std::vector<unsigned> evnm;
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

#if 0
	int j, k, oct, comb_num[8];
	int every;
	double **etas = NULL;

	etas = malloc(8 * sizeof(double));
	for (i = 0; i < 8; i++)
		etas[i] = calloc(ISUB * NUM, sizeof(double));

	for (isub = 0; isub < ISUB; isub++){
		for (irun = 0; irun < NUM; irun++){
			double 
		}
	} /* end of event loop for this event set */

	for (every = 1; every < EVERY_MAX; every++){
		int numevents = NUM*( (ISUB - 1)/every + 1);
		/* more convenient to first calculate then write */
		double corr[8][8];
		for (oct = 0; oct < 8; oct++){
			sprintf (s, "res/out_oct_%i", oct);
			FILE *fout = fopen(s, "a+");
			rsn = 1/sqrt(numevents);
			mean = gsl_stats_mean (etas[oct], every, NUM * ISUB);
			rsdm = rsn * gsl_stats_sd_m (etas[oct], every,
						     NUM * ISUB, mean);
			fprintf (fout, "%.9f\t%.9f\t%.9f\t%.9f\t%.9f\n",
				 rsn, 1.0/numevents,
				 mean, fabs(mean), rsdm);
			fclose(fout);
			/* correlations of the same octant are 1 by definition
			 * so no need to calculate them. set them to
			 * zero so we don't
			 * have problems with plot scaling */
			corr[oct][oct] = 0;
			for (i = oct + 1; i < 8; i++){
				corr[i][oct]
					= corr[oct][i]
					= gsl_stats_correlation (etas[oct],
								 every,
								 etas[i],
								 every,
								 NUM * ISUB);
			}
		}
		sprintf (s, "res/out_cor_every_%i", every);
		FILE *fcor = fopen (s, "a+");
		for (oct = 0; oct < 8; oct++){
			for (i = 0; i < 8; i++){
				fprintf (fcor, "%.9f\t", corr[oct][i]);
			}
			fprintf(fcor, "\n");
		}
		fclose (fcor);
	} /* end of "every" event set loop */
	/* this can be done earlier? */
	for( isub = 0; isub < ISUB; isub++){
		for( irun = 0; irun < NUM; irun++){
			free( all_particles[isub][irun] );
		}
		free(all_particles[isub]);
		free(pcnt[isub]);
		free (psz[isub]);
	}
	for (i = 0; i < 8; i++) free(etas[i]);
	free(etas);
	free(pcnt);
	free (psz);
	free(all_particles);
	time_t t2 = time(NULL);
	printf("finished calculations and wrap-up in %li\n", t2 - t1);

#endif	// comment out main

	return 0;
}
