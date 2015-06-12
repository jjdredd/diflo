#include "dist.hpp"

distribution::distribution(unsigned ny = 100, unsigned nphi = 100)
	: ny(ny), nphi(nphi), f(nphi) {
	d = new unsigned* [ny];
	for(unsigned i = 0; i < ny; i++) d[i] = new unsigned [nphi];
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
				d[ybin][phibin]++;
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
				d[ybin][phibin]++;
			}
		}
	}
}

void distribution::DistTransform(int yn){
	// f.FTrans(d[yn]);
}
