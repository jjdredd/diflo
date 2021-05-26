#ifndef GRID_HPP
#define GRID_HPP

#include <vector>
#include <string>

#include "data.hpp"

struct SymGrid {

	double dims[3];
	double h[3];
	unsigned Nodes[3];

	SymGrid(double, unsigned);
	bool operator==(SymGrid &);
	std::vector<int> space_to_grid(std::vector<double> &);
};

SymGrid MinGrid(const SymGrid &, const SymGrid &);

template <typename T>
struct ScalarGrid {
	ScalarGrid(const SymGrid &);
	~ScalarGrid();
	SymGrid g;
	T ***elem;
};

class ParticleGrid {

public:
	ParticleGrid(SymGrid &, unsigned);
	virtual ~ParticleGrid();

	int Populate(event &);
	void Clear();
	void ShrinkToFit();
	void WriteParticleCount(const std::string &) const;
	bool IsCellValid(unsigned, unsigned, unsigned) const;
	std::vector<particle> & operator() (unsigned, unsigned, unsigned) const;

private:

	SymGrid g;
	int min_particles;
	ScalarGrid< std::vector<particle> > p;
	ScalarGrid<int> p_cnt;

};

// implement grid accumulator class

class ArrayGrid {

public:
	explicit ArrayGrid(const SymGrid &, unsigned);
	~ArrayGrid();
	double & operator() (unsigned, unsigned, unsigned, unsigned);
	unsigned GetCapacity() const;

private:
	SymGrid g;
	int capacity;
	ScalarGrid<double> garray;
};


#endif	// GRID_HPP
