#include "handedness.hpp"

/* calculate eta for event consisting of n particles stored in p
 * write results in etas array (for 8 octants), writes the number
 * of particle triplet combinations that were used to calculate etas */
void EventEta(event &e, double *etas, unsigned *comb_num){
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
		oct = octant(&p[i]);
		for( j = i + 1; j < n; j++){
			if(oct != octant(&p[j])) {
				continue;
			}
			for(k = j + 1; k < n; k++){
				if(oct != octant(&p[k])){
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
