
#include "data.hpp"
#include "grid.hpp"


//
// class Helicity
//

Helicity::Helicity(const SymGrid &g)
	: g(g) {

	helicity = new double**[g.Nodes[0]];
	velocity = new double**[g.Nodes[0]];
	for (unsigned i = 0; i < g.Nodes[0]; i++) {
		helicity[i] = new double*[g.Nodes[1]];
		velocity[i] = new double*[g.Nodes[1]];
		for (unsigned j = 0; j < g.Nodes[1]; j++) {
			helicity[i][j] = new double[g.Nodes[2]]();
			velocity[i][j] = new double[g.Nodes[2]]();
		}
	}
}


Helicity::~Helicity() {
	for (unsigned i = 0; i < g.Nodes[0]; i++) {
		for (unsigned j = 0; j < g.Nodes[1]; j++) {
			delete[] helicity[i][j];
			delete[] velocity[i][j];
		}
		delete[] helicity[i];
		delete[] velocity[i];
	}
	delete[] helicity;
	delete[] velocity;
}


std::vector<double> Helicity::cell_velocity(const std::vector<particle> &pv) {
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


void Helicity::velocity(const ParticleGrid &pg) {

	for (unsigned i = 0; i < g.Nodes[0]; i++) {
		for (unsigned j = 0; j < g.Nodes[1]; j++) {
			for (unsigned k = 0; k < g.Nodes[2]; k++) {
				std::vector<double> v = cell_velocity(pg(i, j, k));
				velocity[i][j][k] = 
			}
		}
	}
}
