#include "dist.hpp"
#include <cmath>

distribution::distribution() : ny(ny), nphi(nphi), f(nphi) {
	d = new unsigned* [ny];
	for(int i = 0; i < ny; i++) d[i] = new unsigned [nphi];
	return;
}

distribution::~distribution(){
	for(int i = 0; i < ny; i++) delete[] d[i];
	delete[] d;
	return;
}

void distribution::DataDist(data& D){
	for(int isub = 0; isub < D.ISUBS; isub++){
		for(int irun = 0; irun < D.NUM; irun++){
			for(unsigned i = 0; i < D.P[isub][irun].size(); i++){
				int ybin = (int) floor(D.P[isub][irun][i].y()
						       * ny) + ny/2;
				int phibin = (int)
					floor(D.P[isub][irun][i].aangle()*nphi)
					+ nphi/2;
				d[ybin][phibin]++;
			}
		}
	}
		
}

void distribution::DataDist(data& D, int type, int charge, bool meson){
	for(int isub = 0; isub < D.ISUBS; isub++){
		for(int irun = 0; irun < D.NUM; irun++){
			for(unsigned i = 0; i < D.P[isub][irun].size(); i++){
				if(!D.P[isub][irun][i].of_type(type,
							       charge, meson))
					continue;

				int ybin = (int) floor(D.P[isub][irun][i].y()
						       * ny) + ny/2;
				int phibin = (int)
					floor(D.P[isub][irun][i].aangle()*nphi)
					+ nphi/2;
				d[ybin][phibin]++;
			}
		}
	}
		
}

void distribution::DistTransform(int yn){
	f.FTrans(d[yn]);
}
