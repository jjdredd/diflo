#ifndef PARTICLE_HPP
#define PARTICLE_HPP
#include <cmath>

#define THRES 0.0001		// comparison threshold
#define MPROT 0.938		// proton mass

struct particle {

	float x, y, z, Px, Py, Pz, P0, b;
	int type, charge;
	bool meson;

	particle();

	particle(float x, float y, float z, float Px,
		 float Py, float Pz, float P0, float b,
		 int type, int charge);

	~particle();

	bool operator==(float *p);
	float y();
	float aangle();
	float p();
	bool of_type(int t, int c, bool m);
};

#endif	// PARTICLE_HPP
