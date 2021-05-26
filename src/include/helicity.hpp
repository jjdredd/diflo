#ifndef HELICITY_HPP
#define HELICITY_HPP

#include "data.hpp"
#include "grid.hpp"


class Helicity {

public:
	explicit Helicity(const SymGrid &);
	virtual ~Helicity();
	void Compute(const ParticleGrid &);
	void Clear();
	bool CopyArray(ArrayGrid &, unsigned) const;
	void WriteOutHelicity(const std::string &) const;

private:
	void ComputeVelocity(const ParticleGrid &);
	void ComputeHelicity();
	std::vector<double> cell_velocity(const std::vector<particle> &) const;
	std::vector<double> cell_rot(unsigned, unsigned, unsigned);
	double cell_helicity(unsigned, unsigned, unsigned);

	ScalarGrid<double> helicity;
	ArrayGrid v;		// velocity
	SymGrid g;
};


#endif	// HELICITY_HPP
