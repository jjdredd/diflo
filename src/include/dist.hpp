#ifndef DIST_HPP
#define DIST_HPP

#include "data.hpp"
#include "fft.hpp"

// distribution over azimuthal angle and rapidity of particles
// of type 'type' and charge 'charge'

class distribution {

public:
	distribution(unsigned ny = 100, unsigned nphi = 100);
	~distribution();
	void DataDist(data&, int type, int charge, bool meson);
	void DataDist(data&);
	void DistTransform();
	// GetType();
	// GetCharge();
	// GetEvent();

private:
	int type, charge;
	unsigned ny, unsigned nphi;
	bool typed;
	unsigned **d;
	fft f;
};

#endif	// DIST_HPP
