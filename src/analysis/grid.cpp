#include <algorithm>
#include <cmath>

#include "grid.hpp"


//
// struct SymGrid
// 

SymGrid::SymGrid(double size, unsigned Nn) {
	// create symmetric grid by size and Nn nodes
	for (unsigned i = 0; i < 3; i++){
		dims[i] = size;
		Nodes[i] = 2*Nn;
		h[i] = 2*dims[i] / Nodes[i];
	}
}

bool SymGrid::operator==(SymGrid &o) {
	for (unsigned i = 0; i < 3; i++){
		if (fabs(dims[i] - o.dims[i]) > h[i]* h[i]){
			return false;
		}
		if (Nodes[i] != o.Nodes[i]) {
			return false;
		}
	}
	return true;	
}

std::vector<int> SymGrid::space_to_grid(std::vector<double> &x) {
	std::vector<int> v(3);
	for (unsigned i = 0; i < 3; i++) {
		double sh_x = x[i] + dims[i];
		v[i] = static_cast<int> (floor(sh_x / h[i]));
		if (v[i] < 0) v[i] = 0;
	}
	return v;
}


// MinGrid is broken do not use
SymGrid MinGrid(const SymGrid &g_1, const SymGrid &g_2) {
	SymGrid g(0, 0);
	// are grids equal?
	// change g.h as well
	for (unsigned i = 0; i < 3; i++){
		g.dims[i] = std::max(g_1.dims[i], g_2.dims[i]);
		g.Nodes[i] = std::max(g_1.Nodes[i], g_2.Nodes[i]);
		g.h[i] = 2*g.dims[i] / g.Nodes[i];
	}
	return g;
}


//
// class ParticleGrid
// 

ParticleGrid::ParticleGrid(SymGrid &g, unsigned mp) : g(g), min_particles(mp) {
	p = new std::vector<particle>**[g.Nodes[0]];
	p_cnt = new int**[g.Nodes[0]];
	for (unsigned i = 0; i < g.Nodes[0]; i++) {
		p[i] = new std::vector<particle>*[g.Nodes[1]];
		p_cnt[i] = new int*[g.Nodes[1]];
		for (unsigned j = 0; j < g.Nodes[1]; j++) {
			p[i][j] = new std::vector<particle>[g.Nodes[2]];
			p_cnt[i][j] = new int[g.Nodes[2]]();
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

void ParticleGrid::Populate(event &e) {
	for (unsigned n = 0; n < e.particles.size(); n++) {
		particle ptcl = e.particles[n];
		std::vector<double> x{ptcl.x, ptcl.y, ptcl.z};
		std::vector<int> v = g.space_to_grid(x);

		// XXX TODO check if coordinates are on grid (oob)
		bool oob_detected = false;
		for (unsigned i = 0; i < 3; i++) {
			if(oob_detected = (v[i] < 0 || v[i] >= g.Nodes[i])) {
				break;
			}
		}
		if (oob_detected) continue;

		p[ v[0] ][ v[1] ][ v[2] ].push_back(ptcl);
		p_cnt[ v[0] ][ v[1] ][ v[2] ] = p[ v[0] ][ v[1] ][ v[2] ].size();
	}
}

void ParticleGrid::Clear() {
	for (unsigned i = 0; i < g.Nodes[0]; i++) {
		for (unsigned j = 0; j < g.Nodes[1]; j++) {
			for (unsigned k = 0; k < g.Nodes[2]; k++) {
				p[i][j][k].clear();
			}
		}
	}
}

void ParticleGrid::ShrinkToFit() {
	for (unsigned i = 0; i < g.Nodes[0]; i++) {
		for (unsigned j = 0; j < g.Nodes[1]; j++) {
			for (unsigned k = 0; k < g.Nodes[2]; k++) {
				p[i][j][k].shrink_to_fit();
			}
		}
	}
}

void ParticleGrid::WriteParticleCount(const std::string &base_path) const {
	for (unsigned j = 0; j < g.Nodes[1]; j++) {
		std::ofstream out_file;
		std::string file_path;
		file_path = base_path + std::to_string(j)
			+ std::string(".txt");
		out_file.open(file_path, std::ofstream::out);
		for (unsigned k = 0; k < g.Nodes[2]; k++) {
			for (unsigned i = 0; i < g.Nodes[0]; i++) {
				out_file << p_cnt[i][j][k] << '\t';
			}
			out_file << std::endl;
		}
		out_file.close();
	}
}
