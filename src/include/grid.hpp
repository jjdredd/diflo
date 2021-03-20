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

	void Populate(event &);
	void Clear();
	void ShrinkToFit();
	void WriteParticleCount(const std::string &) const;

private:
	bool cell_valid(unsigned, unsigned, unsigned);

	SymGrid g;
	int min_particles;
	std::vector<particle> ***p;
	int ***p_cnt;

};


#endif	// GRID_HPP
