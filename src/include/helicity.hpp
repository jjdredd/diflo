#ifndef HELICITY_HPP
#define HELICITY_HPP

#include "data.hpp"
#include "grid.hpp"


class Helicity {

public:
	explicit Helicity(const SymGrid &);
	virtual ~Helicity();
	void Compute(const ParticleGrid &);

private:
	void velocity(const ParticleGrid &);
	void helicity();
	std::vector<double> cell_velocity(const std::vector<particle> &) const;
	double cell_helicity(unsigned, unsigned, unsigned) const;

	double ***helicity;
	ArrayGrid velocity;
	SymGrid g;
};


#endif	// HELICITY_HPP
