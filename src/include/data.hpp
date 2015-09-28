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

inline double dotprod(const particle &, const particle &);
double dist(const particle &, const particle &);
double mixprod(const particle &, const particle &, const particle &);
bool pcompare(const particle, const particle);


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

enum DataVersion {
	HSD_VER_ORIG,
	HSD_VER_COORD,
	HSD_VER_PHSD,
};

class data {

public:
	unsigned ISUBS, NUM, A, NParticles;
	double Elab, Time;
	event **P;
	data(std::ifstream &s, DataVersion v, bool, int);
	~data();
	unsigned NumberOfParticles();
	void readin_particles(std::ifstream &s, bool mesons);

private:
	DataVersion hsd_ver;
	int type, charge;
	bool pick;
	bool parse_input_line(char *str, int *isub, int *irun, particle *p);

};



//
// ALICE event fetcher
//


class ALICEData {

public:
	ALICEData(const char *);
	~ALICEData();
	bool FetchEvent(event &);
	bool FetchNumEvent(event &, unsigned, unsigned);

private:
	unsigned cur_npart, cur_nev;
	std::ifstream s;
	bool parse_hline(char *);
	bool parse_line(char *, particle &);
};

#if 0

//
// another class to represent data simple input
//


class DataSIn {

public:
	unsigned NParticles;
	std::vector<event> Events;
	DataSIn();
	~DataSIn();
	bool FetchEvent(std::ifstream &s, event &e);
	void readin_data(std::ifstream &s);

private:
	DataVersion dat_ver;
	int parse_input_line(char *str, unsigned &e_num, particle *p);
	unsigned EventNum, isub, irun;
	particle p;		// to store particle from the next event

};

#endif	// comment out DataSIn class

#endif	// DATA_HPP
