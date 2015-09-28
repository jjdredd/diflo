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
#define RPAngle 0.0		// angle of the reaction plane

#define SubVolume(p) diant(p)	// choose a volume partition function

unsigned octant(particle &);
unsigned diant(particle &);


void EventEta(event &, std::vector<double> &, std::vector<unsigned> &);


#endif	// HANDED_HPP
