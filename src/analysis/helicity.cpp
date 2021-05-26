
#include "data.hpp"
#include "grid.hpp"
#include "helicity.hpp"


//
// class Helicity
//

Helicity::Helicity(const SymGrid &g)
	: helicity(g), v(g, 3), g(g) {

}


Helicity::~Helicity() {}


std::vector<double> Helicity::cell_velocity(const std::vector<particle> &pv) const {
	double e = 0;
	std::vector<double> vcell(3, 0.0);

	for (const particle &p : pv) {
		vcell[0] += p.Px;
		vcell[1] += p.Py;
		vcell[2] += p.Pz;
		e += p.P0;
	}
	if (pv.size() > 0) {
		for (unsigned i = 0; i < 3; i++) vcell[i] /= e;
	} else {
		for (unsigned i = 0; i < 3; i++) vcell[i] = 0;
	}

	return vcell;
}


void Helicity::ComputeVelocity(const ParticleGrid &pg) {
	for (unsigned i = 0; i < g.Nodes[0]; i++) {
		for (unsigned j = 0; j < g.Nodes[1]; j++) {
			for (unsigned k = 0; k < g.Nodes[2]; k++) {
				std::vector<double> vcell(3, 0.0);
				if (pg.IsCellValid(i, j, k)) {
					vcell = cell_velocity(pg(i, j, k));
				}
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
	for (unsigned i = 1; i < g.Nodes[0]; i++) {
		for (unsigned j = 1; j < g.Nodes[1]; j++) {
			for (unsigned k = 1; k < g.Nodes[2]; k++) {
				helicity.elem[i][j][k] = cell_helicity(i, j, k);
			}
		}
	}
}

void Helicity::Compute(const ParticleGrid &pg) {
	ComputeVelocity(pg);
	ComputeHelicity();
}

bool Helicity::CopyArray(ArrayGrid &out, unsigned n) const {
	if (n >= out.GetCapacity()) return false;
	for (unsigned i = 0; i < g.Nodes[0]; i++) {
		for (unsigned j = 0; j < g.Nodes[1]; j++) {
			for (unsigned k = 0; k < g.Nodes[2]; k++) {
				out(i, j, k, n) = helicity.elem[i][j][k];
			}
		}
	}
	return true;
}


void Helicity::Clear() {
	for (unsigned i = 0; i < g.Nodes[0]; i++) {
		for (unsigned j = 0; j < g.Nodes[1]; j++) {
			for (unsigned k = 0; k < g.Nodes[2]; k++) {
				helicity.elem[i][j][k] = 0;
				for (unsigned n = 0; n < 3; n++) {
					v(i, j, k, n) = 0;
				}
			}
		}
	}
}

void Helicity::WriteOutHelicity(const std::string &base_path) const {
	for (unsigned j = 0; j < g.Nodes[1]; j++) {
		std::ofstream out_file;
		std::string file_path;
		file_path = base_path + std::to_string(j)
			+ std::string(".txt");
		out_file.open(file_path, std::ofstream::out);
		for (unsigned k = 0; k < g.Nodes[2]; k++) {
			for (unsigned i = 0; i < g.Nodes[0]; i++) {
				out_file << helicity.elem[i][j][k] << '\t';
			}
			out_file << std::endl;
		}
		out_file.close();
	}
}
