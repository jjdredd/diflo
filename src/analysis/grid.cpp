#include <algorithm>
#include <cmath>

#include "grid.hpp"


//
// struct Grid
// 

Grid::Grid(double size, unsigned Nn) {
	// create symmetric grid by size and Nn nodes
	for (unsigned i = 0; i < 3; i++){
		dims[i][0] = -size;
		dims[i][1] = size;
		Nodes[i] = 2*Nn;
		h[i] = (dims[i][1] - dims[i][0]) / Nodes[i];
	}
}

bool Grid::operator==(Grid &o) {
	for (unsigned i = 0; i < 3; i++){
		if (fabs(dims[i][0] - o.dims[i][0]) > h[i]* h[i]
		    || fabs(dims[i][1] - o.dims[i][1]) > h[i]* h[i]){
			return false;
		}
	}
	return true;	
}


Grid MinGrid(const Grid &g_1, const Grid &g_2) {
	Grid g;
	for (unsigned i = 0; i < 3; i++){
		g.dims[i][0] = std::max(g_1.dims[i][0], g_2.dims[i][0]);
		g.dims[i][1] = std::min(g_1.dims[i][1], g_2.dims[i][1]);

		g.Nodes[i][0] = std::max(g_1.Nodes[i][0], g_2.Nodes[i][0]);
		g.Nodes[i][1] = std::min(g_1.Nodes[i][1], g_2.Nodes[i][1]);
	}
	return g;
}


//
// class ParticleGrid
// 

ParticleGrid::ParticleGrid(Grid &g) : g(g) {
	p = new void**[g.Nodes[0]];
	valid = new void**[g.Nodes[0]];
	for (unsigned i = 0; i < g.Nodes[1]; i++) {
		p[i] = new void*[g.Nodes[1]];
	}
}

bool ParticleGrid::cell_valid(unsigned i, unsigned j, unsigned k) {


}

void ParticleGrid::Populate(event &e) {

}
