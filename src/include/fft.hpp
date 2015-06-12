#ifndef FFT_HPP
#define FFT_HPP

#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>

// fast Fourier transform using GSL

class fft {

public:
	fft(int n);
	~fft();
	void FTrans(unsigned *a);
	double GetV(unsigned i);
	double GetA(unsigned i);

private:
	// mostly the stuff needed for gsl fft goes here
	gsl_fft_real_wavetable *rwt;
	gsl_fft_halfcomplex_wavetable *hcwt;
	gsl_fft_real_workspace *work;
	unsigned *rda;		// array of real data that is passed in for
				// transform
	unsigned size;		// size of input data
};

#endif	// FFT_HPP
