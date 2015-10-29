#include "momentum.hpp"

// if zdir == true then pick pz > 0 else pick pz < 0
double MeanPx(DataHSD& D, bool zdir){

	double mpx = 0;
	unsigned n = 0;

	for(unsigned isub = 0; isub < D.ISUBS; isub++){
		for(unsigned irun = 0; irun < D.NUM; irun++){
			for(unsigned i = 0;
			    i < D.P[isub][irun].particles.size(); i++){

				if(zdir
				   && (D.P[isub][irun].particles[i].Pz > 0)){

					mpx += D.P[isub][irun].particles[i].Px;
					n++;

				} else if(!zdir
					  && (D.P[isub][irun].particles[i].Pz
					      < 0)){

					mpx += D.P[isub][irun].particles[i].Px;
					n++;

				}

			}
		}
	}
	return mpx/n;
}
