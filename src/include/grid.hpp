#ifndef GRID_HPP
#define GRID_HPP

#include <algorithm>
#include <cmath>

#include "data.hpp"

struct Grid {

	double dims[3][2];	// 0 - start, 1 - end
	double h[3];
	unsigned Nodes[3];

	Grid(double size, unsigned Nn) {
		// create symmetric grid by size and Nn nodes
		for (unsigned i = 0; i < 3; i++){
			dims[i][0] = -size;
			dims[i][1] = size;
			Nodes[i] = Nn;
			h[i] = (dims[i][1] - dims[i][0]) / Nodes[i];
		}
	}

	bool operator==(Grid &o) {
		for (unsigned i = 0; i < 3; i++){
			if (fabs(dims[i][0] - o.dims[i][0]) > h[i]* h[i]
			    || fabs(dims[i][1] - o.dims[i][1]) > h[i]* h[i]){
				return false;
			}
		}
		return true;	
	}
};


Grid MinGrid(const Grid &g_1, const Grid &g_2) {
	Grid g;
	for (unsigned i = 0; i < 3; i++){
		g.dims[i][0] = std::max(g_1.dims[i][0], g_2.dims[i][0]);
		g.dims[i][1] = std::min(g_1.dims[i][1], g_2.dims[i][1]);
	}
	return g;
}


class ParticleGrid {

public:
	ParticleGrid(Grid &);
	virtual ~ParticleGrid();

	void Populate(event &);

private:
	Grid g;
	particle ****p;

};


#endif	// GRID_HPP
