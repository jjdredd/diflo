#include "event.hpp"

event::event(){}

event::~event(){}

void event::add_particle(particle &p){
	particles.push_back(p);
	return;
}

unsigned int event::particle_count(){
	return particles.size();
}
