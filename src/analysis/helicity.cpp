
#include "data.hpp"
#include "grid.hpp"
#include "helicity.hpp"


//
// class Helicity
//

Helicity::Helicity(const SymGrid &g, unsigned N)
	: v(g, N), g(g) {

	helicity = new double**[g.Nodes[0]];
	for (unsigned i = 0; i < g.Nodes[0]; i++) {
		helicity[i] = new double*[g.Nodes[1]];
		for (unsigned j = 0; j < g.Nodes[1]; j++) {
			helicity[i][j] = new double[g.Nodes[2]]();
		}
	}
}


Helicity::~Helicity() {
	for (unsigned i = 0; i < g.Nodes[0]; i++) {
		for (unsigned j = 0; j < g.Nodes[1]; j++) {
			delete[] helicity[i][j];
		}
		delete[] helicity[i];
	}
	delete[] helicity;
}


std::vector<double> Helicity::cell_velocity(const std::vector<particle> &pv) const {
	double e = 0;
	std::vector<double> v(3, 0.0);

	for (const particle &p : pv) {
		v[0] += p.Px;
		v[1] += p.Py;
		v[2] += p.Pz;
		e += p.P0;
	}
	for (unsigned i = 0; i < 3; i++) v[i] /= e;

	return v;
}


void Helicity::ComputeVelocity(const ParticleGrid &pg) {
	for (unsigned i = 0; i < g.Nodes[0]; i++) {
		for (unsigned j = 0; j < g.Nodes[1]; j++) {
			for (unsigned k = 0; k < g.Nodes[2]; k++) {
				std::vector<double> vcell = cell_velocity(pg(i, j, k));
				for (unsigned n = 0; n < 3; n++) {
					v(i, j, k, n) = vcell[n];
				}
			}
		}
	}
}

std::vector<double> Helicity::cell_rot(unsigned i, unsigned j, unsigned k) {
	std::vector<double> rot(3, 0.0);
	rot[0] = (v(i, j, k, 2) - v(i, j - 1, k, 2)) / g.h[1]
		- (v(i, j, k, 1) - v(i, j, k - 1, 1)) / g.h[2];

	rot[1] = (v(i, j, k, 0) - v(i, j, k - 1, 0)) / g.h[2]
		- (v(i, j, k, 2) - v(i - 1, j, k, 2)) / g.h[0];

	rot[2] = (v(i, j, k, 1) - v(i - 1, j, k, 1))/g.h[0]
		- (v(i, j, k, 0) - v(i, j - 1, k, 0))/g.h[1];

	return rot;
}

double Helicity::cell_helicity(unsigned i, unsigned j, unsigned k) {
	// check oob (?)
	std::vector<double> rot = cell_rot(i, j, k);
	double h = 0;
	for (unsigned n = 0; n < 3; n++) {
		h += rot[n] * v(i, j, k, n);
	}
	return h;
}

void Helicity::ComputeHelicity() {
	for (unsigned i = 0; i < g.Nodes[0]; i++) {
		for (unsigned j = 0; j < g.Nodes[1]; j++) {
			for (unsigned k = 0; k < g.Nodes[2]; k++) {
				helicity[i][j][k] = cell_helicity(i, j, k);
			}
		}
	}
}
