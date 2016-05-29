#include "handedness.hpp"
#include <gsl/gsl_statistics_double.h>


#if 0
static bool momemtum_test(particle &a, particle &b, particle &c){
	return (MOMDIST && ((dist(arp[2], arp[1]) > MAX_DIST)
			    || (dist(arp[2], arp[0]) < MIN_DIST)));
}
#endif


//
// Handedness base class
//

Handedness::Handedness() {}

Handedness::~Handedness() {}

// determine the octant in momentum space
unsigned Handedness::sub_volume(particle &p){
	if (p.Px > 0){
		if (p.Py > 0){
			if (p.Pz > 0) return 0;
			else return 1;
		}
		else{
			if (p.Pz > 0) return 2;
			else return 3;
		}
	}
	else{
		if (p.Py > 0){
			if (p.Pz > 0) return 4;
			else return 5;
		}
		else{
			if (p.Pz > 0) return 6;
			else return 7;
		}
	} 
}

/* calculate eta for event consisting of n particles stored in p
 * write results in etas array (for 8 octants), writes the number
 * of particle triplet combinations that were used to calculate etas */
void Handedness::EventEta(event &e, std::vector<double> &etas,
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
		oct = sub_volume(p[i]);
		for( j = i + 1; j < n; j++){
			if(oct != sub_volume(p[j])) continue;
			for(k = j + 1; k < n; k++){
				if(oct != sub_volume(p[k])) continue;
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


//
// HandednessExp
//

HandednessExp::HandednessExp() : RPAngle(0) {}

unsigned HandednessExp::sub_volume(particle &p) {
	if ((-p.Px*tan(RPAngle) + p.Py) > 0) return 0;
	else return 1;
}

double MaxHandedRatio(event& e, double& angle) {
	double prev_ratio = 0;
	HandednessExp H;
	std::vector<double> evet;
	std::vector<unsigned> evnm;
	for (double rpa = 0; rpa < M_PI_2; rpa += 0.05) {
		H.RPAngle = rpa;
		H.EventEta(e, evet, evnm);

		// if ((fabs(evet[0]) < ETA_THRES)
		//     && (fabs(evet[1]) < ETA_THRES)) continue;

		// double ratio = fabs(evet[0] - evet[1])
		// 	/ (fabs(evet[0]) + fabs(evet[1]));

		double ratio =fabs(evet[0])
			+ fabs(evet[1]);

		if (ratio > prev_ratio) {
			prev_ratio = ratio;
			angle = rpa;
		}
	}
	return prev_ratio;
}

void HandedStatExp(std::vector<event>& Events, double angle,
		   std::vector<double>& mean, std::vector<double>& rsdm) {

	std::vector<double> hs[2];
	HandednessExp H;

	mean.reserve(2);
	rsdm.reserve(2);
	for (auto &e : Events) {
		std::vector<double> eta;
		std::vector<unsigned> num;
		H.RPAngle = e.RPA + angle;
		H.EventEta(e, eta, num);
		hs[0].push_back(eta[0]);
		hs[1].push_back(eta[1]);
	}
	for (unsigned i = 0; i < 2; i++) {
		double rsn = 1 / sqrt(Events.size());
		mean[i] = gsl_stats_mean (&hs[i][0], 1, Events.size());
		rsdm[i] = rsn * gsl_stats_sd_m (&hs[i][0], 1,
						Events.size(), mean[i]);
	}
}
