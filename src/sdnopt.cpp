/*****************************************************************************
 * This program is free software: you can redistribute it and/or modify	     *
 * it under the terms of the GNU General Public License version 2.	     *
 * 									     *
 * Author: Rahim Usubov <usubov@theor.jinr.ru>				     *
 * 									     *
 * GNU GSL is required to build this program. Build detailes are provided in *
 * the Makefile.							     *
 *****************************************************************************/

#include <iostream>
#include <cstdio>

#include <gsl/gsl_statistics_double.h>
/* HSD and PHSD seem to have slightly different output formats,
 * to take this into account, a PHSD preprocessor variable is used */
/* #define PHSD */
#define MIN_DIST 1.0
#define MAX_DIST 1.0
#define MOMDIST 0		/* if zero skip momentum test */
#define norm(a) sqrt(dotprod(a, a))

struct particle {
	double P[4]; 			/* 0 - E, 3 - Pz */
	/* see hsd manual for the meaning and values of these variables */
	int charge, isub, irun, type;
};

inline int is_pion_event(struct particle *p, int isub, int irun){
	return ( (p->type == 1)
		 && (p->isub == isub)
		 && (p->irun == irun) );
}

int quadrant(struct particle *p){
	if( p->P[1] > 0 ){
		if( p->P[2] > 0 ){
			if( p->P[3] > 0 ) return 0;
			else return 1;
		}
		else{
			if( p->P[3] > 0 ) return 2;
			else return 3;
		}
	}
	else{
		if( p->P[2] > 0 ){
			if( p->P[3] > 0 ) return 4;
			else return 5;
		}
		else{
			if( p->P[3] > 0 ) return 6;
			else return 7;
		}
	} 
}

inline double dotprod( struct particle *p1, struct particle *p2){
	double sdp = 0;
	int i;
	for( i = 1; i <=3; i++) sdp += p1->P[i]*p2->P[i]; 
	return sdp;
}

inline double dist( struct particle *p1, struct particle *p2){
	struct particle p;
	int i;
	for( i = 1; i <= 3; i++) p.P[i] = p1->P[i] - p2->P[i];
	return norm(&p);
}

/* accepts sorted by inplace_3sort() array */
double mixprod( struct particle **p){
	struct particle r;
	r.P[1] = p[1]->P[2]*p[0]->P[3] - p[1]->P[3]*p[0]->P[2];
	r.P[2] = p[1]->P[3]*p[0]->P[1] - p[1]->P[1]*p[0]->P[3];
	r.P[3] = p[1]->P[1]*p[0]->P[2] - p[1]->P[2]*p[0]->P[1];
	return dotprod(p[2], &r);
}

int pcompare(struct particle *a, struct particle *b){
	double c = dotprod(a, a) - dotprod(b, b);
	if(c > 0) return 1;
	else{
		if(c < 0) return -1;
		else return 0;
	}
}

/* not used atm */
void inplace_3sort(struct particle **arp){
	int min, max, mid, i;		/* mid - hardest lane */
	struct particle *p[3];
	min = max = mid = 0;
	for( i = 0; i < 3; i++){
		if( norm(arp[i]) < norm(arp[min]) ) min = i;
		if( norm(arp[i]) > norm(arp[max]) ) max = i;
	}
	/* find mid :) */
	for( i = 0; i < 3; i++){
		if( (i != min) && (i != max) ){
			mid = i;
			break;
		}
	}
	/* in-place */
	p[0] = arp[max];
	p[1] = arp[mid];
	p[2] = arp[min];
	/* p[0] > p[1] > p[2] */
	for( i = 0; i < 3; i++) arp[i] = p[i];
	/* arp[0] > arp[1] > arp[2] */
	return;
}

/* calculate eta for event consisting of n particles stored in p
 * write results in etas array (for 8 octants), writes the number
 * of particle triplet combinations that were used to calculate etas */
void EventEta(particle *p, unsigned n, double *etas, unsigned *comb_num){
	struct particle *arp[3];
	double eta[8], abseta[8], mxpr;
	unsigned oct, i, j, k;
	memset(eta, 0, sizeof(eta));
	memset(abseta, 0, sizeof(abseta));
	memset(comb_num, 0, sizeof(comb_num));

	qsort(p, n, sizeof(struct particle), &pcompare);

	/* traverse all particles in this
	 * event (isub, inum) */
	for ( i = 0; i < n; i++){
		oct = quadrant(&p[i]);
		for( j = i + 1; j < n; j++){
			if(oct != quadrant(&p[j])) {
				continue;
			}
			for(k = j + 1; k < n; k++){
				if(oct != quadrant(&p[k])){
					continue;
				}
				/* first sort particles */
				/* always i < j < k
				 * -> if presorted then
				 * p_i < p_j < p_k */
				arp[0] = &p[k];
				arp[1] = &p[j];
				arp[2] = &p[i];

				/* arp[0] > arp[1] > arp[2] */
				/* continue if particles are
				 * at wrong distance */
				if(MOMDIST
				   && ((dist(arp[2], arp[1]) > MAX_DIST)
				       || (dist(arp[2], arp[0]) < MIN_DIST))) {
					continue;
				}

				mxpr = mixprod(arp);
				eta[oct] += mxpr;
				abseta[oct] += fabs(mxpr);
				comb_num[oct]++;
			}
		}
	} /* end of particle loop */
	/* update mean if we had particles
	 * in octant and event */
	for (oct = 0; oct < 8; oct++){
		if (comb_num[oct]) {
			etas[oct] = eta[oct] / abseta[oct];
		}
	}

}

int main(int argc, char** argv){
  
	if (argc < 3){
		printf ("not enough arguments\n"
			"usage:\n"
			"first arg fort.301\n"
			"second arg input\n");
		return -1;
	}
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
	return 0;
}
