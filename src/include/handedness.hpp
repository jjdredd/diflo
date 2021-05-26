#ifndef HANDED_HPP
#define HANDED_HPP

#include "data.hpp"
#include "grid.hpp"
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
		      std::vector<unsigned> &);

private:
	virtual unsigned sub_volume(particle &); // octants
};


//
// HandednessExp class (diants)
//

class HandednessExp : public Handedness {

public:
	HandednessExp();
	double RPAngle;

private:
	virtual unsigned sub_volume(particle &p);
};

double MaxHandedRatio(event&, double&);


//
// Handedness on a spacial grid
// 

class HandednessGrid {

public:
	HandednessGrid(const SymGrid &);
	virtual ~HandednessGrid();
	void Compute(const ParticleGrid &);
	void WriteOutHandedness(const std::string &) const;
	void Clear();
	bool CopyArray(ArrayGrid &, unsigned) const;

private:
	double compute_cell_hand(const std::vector<particle> &) const;

	SymGrid g;
	double ***hand;

};



#endif	// HANDED_HPP
