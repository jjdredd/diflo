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
	float GetV(unsigned i);
	float GetA(unsigned i);

private:
	// mostly the stuff needed for gsl fft goes here
	gsl_fft_real_wavetable *rwt;
	gsl_fft_halfcomplex_wavetable *hcwt;
	gsl_fft_read_workspace *work;
	unsigned *rda;		// array of real data that is passed in for
				// transform
};

#endif	// FFT_HPP
