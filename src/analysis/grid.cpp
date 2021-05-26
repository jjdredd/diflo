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
// ScalarGrid
//

template <typename T>
ScalarGrid::ScalarGrid(const SymGrid &g)
	: g(g) {
	elem = new T**[g.Nodes[0]];
	for (unsigned i = 0; i < g.Nodes[0]; i++) {
		elem[i] = new T*[g.Nodes[1]];
		for (unsigned j = 0; j < g.Nodes[1]; j++) {
			elem[i][j] = new T[g.Nodes[2]];
		}
	}
}

template <typename T>
ScalarGrid::~ScalarGrid() {
	for (unsigned i = 0; i < g.Nodes[0]; i++) {
		for (unsigned j = 0; j < g.Nodes[1]; j++) {
			delete[] elem[i][j];
		}
		delete[] elem[i];
	}
	delete[] elem;
}

template<> struct ScalarGrid<int>;
template<> struct ScalarGrid<double>;
template<> struct ScalarGrid< std::vector<particle> >;


//
// class ParticleGrid
// 

ParticleGrid::ParticleGrid(SymGrid &g, unsigned mp)
	: g(g), min_particles(mp), p(g), p_cnt(g) {
}

ParticleGrid::~ParticleGrid() {}

bool ParticleGrid::IsCellValid(unsigned i, unsigned j, unsigned k) const {

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

int ParticleGrid::Populate(event &e) {
	int oob_count = 0;
	for (unsigned n = 0; n < e.particles.size(); n++) {
		particle ptcl = e.particles[n];
		std::vector<double> x{ptcl.x, ptcl.y, ptcl.z};
		std::vector<int> v = g.space_to_grid(x);

		// XXX TODO check if coordinates are on grid (oob)
		bool oob_detected = false;
		for (unsigned i = 0; i < 3; i++) {
			if( (oob_detected = (v[i] < 0 || v[i] >= g.Nodes[i])) ) {
				break;
			}
		}
		if (oob_detected) {
			oob_count++;
			continue;
		}

		p[ v[0] ][ v[1] ][ v[2] ].push_back(ptcl);
		p_cnt[ v[0] ][ v[1] ][ v[2] ] = p[ v[0] ][ v[1] ][ v[2] ].size();
	}
	return oob_count;
}

void ParticleGrid::Clear() {
	for (unsigned i = 0; i < g.Nodes[0]; i++) {
		for (unsigned j = 0; j < g.Nodes[1]; j++) {
			for (unsigned k = 0; k < g.Nodes[2]; k++) {
				p[i][j][k].clear();
				p_cnt[i][j][k] = 0;
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


std::vector<particle> & ParticleGrid::operator() (unsigned i, unsigned j, unsigned k) const {
	return p[i][j][k];
}



//
// class ArrayGrid
//

ArrayGrid::ArrayGrid(const SymGrid &g, unsigned size) : g(g), capacity(size) {
	garray = new double***[g.Nodes[0]];
	for (unsigned i = 0; i < g.Nodes[0]; i++) {
		garray[i] = new double**[g.Nodes[1]];
		for (unsigned j = 0; j < g.Nodes[1]; j++) {
			garray[i][j] = new double*[g.Nodes[2]];
			for (unsigned k = 0; k < g.Nodes[2]; k++) {
				garray[i][j][k] = new double[capacity];
			}
		}
	}
}


ArrayGrid::~ArrayGrid() {
	for (unsigned i = 0; i < g.Nodes[0]; i++) {
		for (unsigned j = 0; j < g.Nodes[1]; j++) {
			for (unsigned k = 0; k < g.Nodes[2]; k++) {
				delete[] garray[i][j][k];
			}
			delete[] garray[i][j];
		}
		delete[] garray[i];
	}
	delete[] garray;
}

double & ArrayGrid::operator() (unsigned i, unsigned j, unsigned k, unsigned n) {
	if (n >= capacity) {
		std::cerr << "OOB in ArrayGrid::operator()" << std::endl;
	}
	return garray[i][j][k][n];
}

unsigned ArrayGrid::GetCapacity() const {
	return capacity;
}
