#ifndef HELICITY_HPP
#define HELICITY_HPP

#include "data.hpp"
#include "grid.hpp"


// make single-responsibility classes!
class Helicity {

public:
	Helicity();
	virtual ~Helicity();
	void Compute(const event &);

private:
	void check_bounds(const event &);
	void helicity(const event &);
	void handedness(const event &);

	double ***helicity;
	Grid g;
	
};


#endif	// HELICITY_HPP
