#include "dist.hpp"

distribution::distribution(unsigned ny, unsigned nphi)
	: ny(ny), nphi(nphi), f(nphi) {
	d = new double* [ny];
	for(unsigned i = 0; i < ny; i++) d[i] = new double [nphi];
	return;
}

distribution::~distribution(){
	for(unsigned i = 0; i < ny; i++) delete[] d[i];
	delete[] d;
	return;
}

void distribution::DataDist(data& D){

	double rapidity, angle;
	int ybin, phibin;

	for(int isub = 0; isub < D.ISUBS; isub++){
		for(int irun = 0; irun < D.NUM; irun++){
			for(unsigned i = 0;
			    i < D.P[isub][irun].particles.size(); i++){

				rapidity = D.P[isub][irun].particles[i].rapid();
				angle = D.P[isub][irun].particles[i].aangle();
				ybin = (int) floor(rapidity * ny/2) + ny/2;
				phibin = (int) floor(angle * nphi/M_PI)
					+ nphi/2;
				d[ybin][phibin] += 1.0;
			}
		}
	}
}

void distribution::DataDist(data& D, int type, int charge, bool meson){

	double rapidity, angle;
	int ybin, phibin;

	for(int isub = 0; isub < D.ISUBS; isub++){
		for(int irun = 0; irun < D.NUM; irun++){
			for(unsigned i = 0;
			    i < D.P[isub][irun].particles.size(); i++){

				if(!D.P[isub][irun].particles[i]
				   .of_type(type, charge, meson))
					continue;

				rapidity = D.P[isub][irun].particles[i].rapid();
				angle = D.P[isub][irun].particles[i].aangle();
				ybin = (int) floor(rapidity * ny/2) + ny/2;
				phibin = (int) floor(angle * nphi/M_PI)
					+ nphi/2;
				d[ybin][phibin] += 1.0;
			}
		}
	}
}

void distribution::DistTransform(double rapidity){
	int ybin = (int) floor(rapidity * ny/2) + ny/2;
	f.FTrans(d[ybin]);
}

void distribution::PrintFlows(std::ostream &s){
	s << "v_2 = " << f.GetV(2)
	  << "\na_2 = " << f.GetA(2) << '\n';
}
