#ifndef HANDED_HPP
#define HANDED_HPP

#include "data.hpp"
#include <cmath>
#include <cstring>
#include <vector>
#include <algorithm>

// move this crap into the class?
#define MIN_DIST 1.0
#define MAX_DIST 1.0
#define MOMDIST 0		/* if zero skip momentum test */


//
// Handedness class
//
 
class Handedness {

public:
	Handedness();
	virtual ~Handedness();
	void EventEta(event &, std::vector<double> &,
			      std::vector<unsigned> &, double);

private:
	virtual unsigned sub_volume(particle &); // octants
};


//
// HandednessExp class (diants)
//

class HandednessExp : public Handedness {

public:
	Handedness();
	~Handedness();

private:
	virtual unsigned sub_volume(particle &p, double angle = 0);
};

#endif	// HANDED_HPP
