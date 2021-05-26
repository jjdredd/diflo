#include "handedness.hpp"


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



//
// Handedness on a spacial grid
// 

HandednessGrid::HandednessGrid(const SymGrid &g) : g(g), hand(g) {}

HandednessGrid::~HandednessGrid() {}

void HandednessGrid::WriteOutHandedness(const std::string &base_path) const {
	for (unsigned j = 0; j < g.Nodes[1]; j++) {
		std::ofstream out_file;
		std::string file_path;
		file_path = base_path + std::to_string(j)
			+ std::string(".txt");
		out_file.open(file_path, std::ofstream::out);
		for (unsigned k = 0; k < g.Nodes[2]; k++) {
			for (unsigned i = 0; i < g.Nodes[0]; i++) {
				out_file << hand.elem[i][j][k] << '\t';
			}
			out_file << std::endl;
		}
		out_file.close();
	}
}

double HandednessGrid::compute_cell_hand(const std::vector<particle> &p) const {
	double eta = 0, abseta = 0;
	std::vector<particle> ptcls(p.size());

	std::copy(p.begin(), p.end(), ptcls.begin());
	std::sort(ptcls.begin(), ptcls.end(), pcompare);
	unsigned size = ptcls.size();

	for (unsigned i = 0; i < size; i++){
		for(unsigned j = i + 1; j < size; j++){
			for(unsigned k = j + 1; k < size; k++){
				double mxpr = mixprod(p[k], p[j], p[i]);
				eta += mxpr;
				abseta += fabs(mxpr);
			}
		}
	}
	return eta / abseta;
}

void HandednessGrid::Compute(const ParticleGrid &pg) {
	for (unsigned i = 0; i < g.Nodes[0]; i++) {
		for (unsigned j = 0; j < g.Nodes[1]; j++) {
			for (unsigned k = 0; k < g.Nodes[2]; k++) {
				if(!pg.IsCellValid(i, j, k)) continue;
				hand.elem[i][j][k] = compute_cell_hand(pg(i, j, k));
			}
		}
	}
}


bool HandednessGrid::CopyArray(ArrayGrid &out, unsigned n) const {
	if (n >= out.GetCapacity()) return false;
	for (unsigned i = 0; i < g.Nodes[0]; i++) {
		for (unsigned j = 0; j < g.Nodes[1]; j++) {
			for (unsigned k = 0; k < g.Nodes[2]; k++) {
				out(i, j, k, n) = hand.elem[i][j][k];
			}
		}
	}
	return true;
}


void HandednessGrid::Clear() {
	for (unsigned i = 0; i < g.Nodes[0]; i++) {
		for (unsigned j = 0; j < g.Nodes[1]; j++) {
			for (unsigned k = 0; k < g.Nodes[2]; k++) {
				hand.elem[i][j][k] = 0;
			}
		}
	}
}
