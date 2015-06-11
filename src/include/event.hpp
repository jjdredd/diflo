#ifndef EVENT_HPP
#define EVENT_HPP
#include <vector>
#include "particle.hpp"

struct event {

	std::vector<particle> particles;
	int A;
	event();
	~event();
	void add_particle(particle &p, bool meson);
};

#endif	// EVENT_HPP
