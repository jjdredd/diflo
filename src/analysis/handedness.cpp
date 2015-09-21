#include "handedness.hpp"


#if 0
static bool momemtum_test(particle &a, particle &b, particle &c){
	return (MOMDIST && ((dist(arp[2], arp[1]) > MAX_DIST)
			    || (dist(arp[2], arp[0]) < MIN_DIST)));
}
#endif

/* calculate eta for event consisting of n particles stored in p
 * write results in etas array (for 8 octants), writes the number
 * of particle triplet combinations that were used to calculate etas */
void EventEta(event &e, std::vector<double> &etas,
	      std::vector<unsigned> &comb_num){

	double eta[8], abseta[8], mxpr;
	unsigned oct, i, j, k;
	std::vector<particle> p = e.particles;
	unsigned n = p.size();
	memset(eta, 0, sizeof(eta));
	memset(abseta, 0, sizeof(abseta));
	comb_num.assign(8, 0);
	etas.assign(8, 0);

	std::sort(p.begin(), p.end(), pcompare);

	/* traverse all particles in this
	 * event (isub, inum) */
	for (i = 0; i < n; i++){
		oct = SubVolume(p[i]);
		for( j = i + 1; j < n; j++){
			if(oct != SubVolume(p[j])) continue;
			for(k = j + 1; k < n; k++){
				if(oct != SubVolume(p[k])) continue;
				/* first sort particles */
				/* always i < j < k
				 * -> if presorted then
				 * p_i < p_j < p_k */

#if 0
				// arp[0] = &p[k];
				// arp[1] = &p[j];
				// arp[2] = &p[i];

				/* arp[0] > arp[1] > arp[2] */
				/* continue if particles are
				 * at wrong distance */
				if (momentum_test(&p[k], &p[j], &p[i]))
					continue;
#endif

				mxpr = mixprod(p[k], p[j], p[i]);
				eta[oct] += mxpr;
				abseta[oct] += fabs(mxpr);
				comb_num[oct]++;
			}
		}
	} /* end of particle loop */
	/* update mean if we had particles
	 * in octant and event */
	for (oct = 0; oct < 8; oct++){
		if (comb_num[oct] && eta[oct]) {
			etas[oct] = eta[oct] / abseta[oct];
		}
	}
}
