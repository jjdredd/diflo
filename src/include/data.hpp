#ifndef DATA_HPP
#define DATA_HPP
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdio>
#include <cmath>
#include "particle.hpp"
#include "event.hpp"


class data {

public:
	int ISUBS, NUM, A;
	float Elab, Time;
	event **P;
	data(std::ifstream &s);
	~data();
	void report_pnum(std::ostream &os);
	void readin_particles(std::ifstream &s, bool mesons);
};

#endif	// DATA_HPP
