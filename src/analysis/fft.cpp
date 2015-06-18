#include "fft.hpp"

fft::fft(int n) : size(n) {
	work = gsl_fft_real_workspace_alloc(size);
	rwt = gsl_fft_real_wavetable_alloc(size);
	rda = new double [size];
	uphc = new double [2*size];
}

fft::~fft(){
	gsl_fft_real_workspace_free(work);
	gsl_fft_real_wavetable_free(rwt);
	delete[] rda;
	delete[] uphc;
}

void fft::FTrans(double *a){

	// make simple accounting here:
	// - make a flag transformed / to be transformed?
	// - errors / return values?

	memcpy(rda, a, size * sizeof(double));
	gsl_fft_real_transform(rda, 1, size, rwt, work);
	gsl_fft_halfcomplex_unpack(rda, uphc, 1, size);
}

double fft::GetV(unsigned i){
	if(2*i < size) return uphc[2*i];
	else return 0;		// error !!!!!
}

double fft::GetA(unsigned i){
	if(2*i + 1 < size) return uphc[2*i + 1];
	else return 0;		// error !!!!!
}
