#include "dist.hpp"

distribution::distribution(unsigned ny, unsigned nphi)
	: ny(ny), nphi(nphi), f(nphi) {
	d = new double* [ny];
	for(unsigned i = 0; i < ny; i++) d[i] = new double [nphi]();
	return;
}

distribution::~distribution(){
	for(unsigned i = 0; i < ny; i++) delete[] d[i];
	delete[] d;
	return;
}

unsigned distribution::BinY(double rapidity) {
	int ybin = (int) floor(rapidity * ny/2) + ny/2;
	if (ybin == ny) ybin = ny -1;
	if (ybin == -1) ybin = 0;
	return ybin;
}

unsigned distribution::BinPhi(double angle) {

	int phibin = (int) floor(angle * nphi/M_PI) + nphi/2;
	// this can happen
	if (phibin == nphi) phibin = nphi - 1;
	if (phibin == -1) phibin = 0;
	return phibin;
}

void distribution::DataDist(DataHSD& D){

	double rapidity, angle;
	unsigned ybin, phibin, pnum = D.NumberOfParticles();

	for(unsigned isub = 0; isub < D.ISUBS; isub++){
		for(unsigned irun = 0; irun < D.NUM; irun++){
			for(unsigned i = 0;
			    i < D.P[isub][irun].particles.size(); i++){

				rapidity = D.P[isub][irun].particles[i].rapid();
				angle = D.P[isub][irun].particles[i].aangle();

				ybin = BinY(rapidity);
				phibin = BinPhi(angle);

				d[ybin][phibin] += 1.0;
			}
		}
	}

	// now normalize to get a distribution
	for (unsigned i = 0; i < ny; i++)
		for (unsigned j = 0; j < nphi; j++)
			d[i][j] /= pnum;
}

void distribution::DataDist(DataHSD& D, int type, int charge, bool meson){

	double rapidity, angle;
	int ybin, phibin, pnum = D.NumberOfParticles();

	for(unsigned isub = 0; isub < D.ISUBS; isub++){
		for(unsigned irun = 0; irun < D.NUM; irun++){
			for(unsigned i = 0;
			    i < D.P[isub][irun].particles.size(); i++){

				if(!D.P[isub][irun].particles[i]
				   .of_type(type, charge, meson))
					continue;

				rapidity = D.P[isub][irun].particles[i].rapid();
				angle = D.P[isub][irun].particles[i].aangle();

				ybin = BinY(rapidity);
				phibin = BinPhi(angle);

				d[ybin][phibin] += 1.0;
			}
		}
	}

	// now normalize to get a distribution
	for (unsigned i = 0; i < ny; i++)
		for (unsigned j = 0; j < nphi; j++)
			d[i][j] /= pnum;
}

void distribution::DistTransform(double rapidity){
	int ybin = BinY(rapidity);
	f.FTrans(d[ybin]);
}

void distribution::PrintFlows(std::ostream &s){
	s << "v_2 = " << f.GetV(2)
	  << "\na_2 = " << f.GetA(2) << '\n';
}

void distribution::WriteDistr(std::ostream &s){
	int ybin;
	const double dphi = M_PI/nphi;
	for (unsigned i = 0; i < nphi; i++){
		for (double r = -1.0; r < 1.0; r += 0.1){
			ybin = BinY(r);
			s << i*dphi - M_PI_2 << '\t' << r << '\t'
			  << d[ybin][i] << std::endl;
		}
		s << std::endl;
	}
}

void distribution::WriteDistrFT(std::ostream &s) {
	for (unsigned i = 0; i < nphi; i++)
		s << i << '\t' << f.GetV(i) << '\t'
		  << f.GetA(i) << std::endl;
}

double RPA_by_multip(event& e) {
	double x = 0, y = 0;
	for (particle p : e.particles) {
		x += p.Px / p.Pt();	// cos()
		y += p.Py / p.Pt();	// sin()
	}
	return atan(y/x);
}

double RPA_by_Pt(event& e) {
	double x = 0, y = 0;
	for (particle p : e.particles) {
		x += p.Px;	// cos()
		y += p.Py;	// sin()
	}
	return atan(y/x);
}
