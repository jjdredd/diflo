#ifndef GRID_HPP
#define GRID_HPP

#include <vector>

#include "data.hpp"

struct Grid {

	double dims[3][2];	// 0 - start, 1 - end
	double h[3];
	unsigned Nodes[3];

	Grid(double, unsigned);
	bool operator==(Grid &);
};

Grid MinGrid(const Grid &, const Grid &);


class ParticleGrid {

public:
	ParticleGrid(Grid &, unsigned);
	virtual ~ParticleGrid();

	void Populate(event &);

private:
	bool cell_valid(unsigned, unsigned, unsigned);
	std::vector<unsigned> space_to_grid(std::vector<unsigned> &);

	Grid g;
	unsigned min_particles;
	std::vector<particle> ***p;
	unsigned char ***p_cnt;

};


#endif	// GRID_HPP
