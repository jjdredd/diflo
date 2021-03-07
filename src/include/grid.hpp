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
	ParticleGrid(Grid &);
	virtual ~ParticleGrid();

	void Populate(event &);

private:
	bool cell_valid(unsigned, unsigned, unsigned);

	Grid g;
	std::vector<particle> ***p;
	unsigned char ***valid;

};


#endif	// GRID_HPP
