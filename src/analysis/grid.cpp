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
		if (Nodes[i] != o.Nodes[i]) {
			return false;
		}
	}
	return true;	
}


// MinGrid is broken do not use
Grid MinGrid(const Grid &g_1, const Grid &g_2) {
	Grid g;
	// are grids equal?
	// change g.h as well
	for (unsigned i = 0; i < 3; i++){
		g.dims[i][0] = std::max(g_1.dims[i][0], g_2.dims[i][0]);
		g.dims[i][1] = std::min(g_1.dims[i][1], g_2.dims[i][1]);

		g.Nodes[i][0] = std::max(g_1.Nodes[i][0], g_2.Nodes[i][0]);
		g.Nodes[i][1] = std::min(g_1.Nodes[i][1], g_2.Nodes[i][1]);

1		h[i] = (dims[i][1] - dims[i][0]) / Nodes[i];
	}
	return g;
}


//
// class ParticleGrid
// 

ParticleGrid::ParticleGrid(Grid &g, unsigned mp) : g(g), min_partiles(mp) {
	p = new void**[g.Nodes[0]];
	p_cnt = new void**[g.Nodes[0]];
	for (unsigned i = 0; i < g.Nodes[0]; i++) {
		p[i] = new void*[g.Nodes[1]];
		p_cnt[i] = new void*[g.Nodes[1]];
		for (unsigned j = 0; j < g.Nodes[1]; j++) {
			p[i][j] = new vector<particle>[g.Nodes[2]];
			p_cnt[i][j] = new unsigned[g.Nodes[2]];
		}
	}
}

ParticleGrid::~ParticleGrid() {
	for (unsigned i = 0; i < g.Nodes[0]; i++) {
		for (unsigned j = 0; j < g.Nodes[1]; j++) {
			delete[] p[i][j];
			delete[] p_cnt[i][j];
		}
		delete[] p[i];
		delete[] p_cnt[i];
	}
	delete[] p;
	delete[] p_cnt;
}

bool ParticleGrid::cell_valid(unsigned i, unsigned j, unsigned k) {

	if (p_cnt[i][j][k] < min_particles) return false;
	for (unsigned n = i; n <= i + 1; n += 2) {
		if (p_cnt[n][j][k] < min_particles) return false;
	}
	for (unsigned n = j; n <= j + 1; n += 2) {
		if (p_cnt[i][n][k] < min_particles) return false;
	}
	for (unsigned n = k; n <= k + 1; n += 2) {
		if (p_cnt[i][j][n] < min_particles) return false;
	}
	return true;
}

std::vector<unsigned> ParticleGrid::space_to_grid(std::vector<unsigned> &x) {
	std::vector<unsigned> v(3);
	for (unsigned i = 0; i < 3; i++) {
		sh_x = x[i] + dims ???
		v[i] = static_cast<unsigned>
	}
}

void ParticleGrid::Populate(event &e) {

}
