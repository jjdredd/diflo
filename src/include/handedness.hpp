#ifndef HANDED_HPP
#define HANDED_HPP

#include "data.hpp"
#include <cmath>
#include <cstring>
#include <vector>
#include <algorithm>

#define MIN_DIST 1.0
#define MAX_DIST 1.0
#define MOMDIST 0		/* if zero skip momentum test */

#define SubVolume(p, a) diant(p, a)	// choose a volume partition function

unsigned octant(particle &);
unsigned diant(particle &, double);


void EventEta(event &, std::vector<double> &,
	      std::vector<unsigned> &, double);


#endif	// HANDED_HPP
