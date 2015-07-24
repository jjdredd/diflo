#ifndef DATA_HPP
#define DATA_HPP
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdio>
#include <cmath>

//
// struct particle
// to represent a single particle
//


#define THRES 0.0001		// comparison threshold
#define MPROT 0.938		// proton mass

struct particle {

	double x, y, z, Px, Py, Pz, P0, b;
	int type, charge;
	bool meson;

	particle();

	particle(double x, double y, double z, double Px,
		 double Py, double Pz, double P0, double b,
		 int type, int charge, bool meson);

	~particle();

	bool operator==(double *p);
	double rapid();
	double aangle();
	double p();
	bool of_type(int t, int c, bool m);
};

//
// struct event
// to represent a collision event
//

struct event {

	std::vector<particle> particles;
	int A;
	event();
	~event();
	void add_particle(particle &p);
	unsigned particle_count();
};

//
// data class
// represents a whole set of input data
//

enum HSDVersion {
	HSD_VER_ORIG,
	HSD_VER_COORD,
	HSD_VER_PHSD
};

class data {

public:
	int ISUBS, NUM, A;
	double Elab, Time, NParticles;
	event **P;
	data(std::ifstream &s, HSDVersion v);
	~data();
	unsigned NumberOfParticles();
	void readin_particles(std::ifstream &s, bool mesons);

private:
	HSDVersion hsd_ver;
	bool parse_input_line(char *str, int *isub, int *irun, particle *p);

};

#endif	// DATA_HPP
