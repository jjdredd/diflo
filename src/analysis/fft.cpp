#include "fft.hpp"

fft::fft(int n) : size(n) {
	work = gsl_fft_real_workspace_alloc(size);
	rwt = gsl_fft_real_wavetable_alloc(size);
	rda = new double [size];
}

fft::~fft(){
	gsl_fft_real_workspace_free(work);
	gsl_fft_real_wavetable_free(rwt);
	delete[] rda;
}

void fft::FTrans(double *a){

	// make simple accounting here:
	// - store results in a seperate array?
	// - make a flag transformed / to be transformed?
	// - errors / return values?

	memcpy(rda, a, size * sizeof(double));
	gsl_fft_real_transform(rda, 1, size, rwt, work);
	// perhaps unpack here?

}

double fft::GetV(unsigned i){}

double fft::GetA(unsigned i){}
