#ifndef GRID_HPP
#define GRID_HPP

#include <vector>
#include <string>

#include "data.hpp"

struct SymGrid {

	double dims[3];
	double h[3];
	unsigned Nodes[3];

	SymGrid(double, unsigned);
	bool operator==(SymGrid &);
	std::vector<int> space_to_grid(std::vector<double> &);
};

SymGrid MinGrid(const SymGrid &, const SymGrid &);


class ParticleGrid {

public:
	ParticleGrid(SymGrid &, unsigned);
	virtual ~ParticleGrid();

	int Populate(event &);
	void Clear();
	void ShrinkToFit();
	void WriteParticleCount(const std::string &) const;
	bool IsCellValid(unsigned, unsigned, unsigned) const;
	std::vector<particle> & operator() (unsigned, unsigned, unsigned) const;

private:

	SymGrid g;
	int min_particles;
	std::vector<particle> ***p;
	int ***p_cnt;

};

// implement grid accumulator class


#endif	// GRID_HPP
