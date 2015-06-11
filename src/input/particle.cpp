#include "particle.hpp"
#include <cmath>

particle::particle(){}

particle::particle(float x, float y, float z, float Px,
		   float Py, float Pz, float P0, float b,
		   int type, int charge, bool meson)
	: x(x), y(y), z(z), Px(Px), Py(Py), Pz(Pz), P0(P0), b(b),
	  type(type), charge(charge), meson(meson){}

particle::~particle(){}

bool particle::operator==(float *p){
	return ((fabsf(this->P0 - p[0]) < THRES)
		&& (fabsf(this->Px - p[1]) < THRES)
		&& (fabsf(this->Py - p[2]) < THRES)
		&& (fabsf(this->Pz - p[3]) < THRES));
}

float particle::y(){
	return atanf(Pz/P0);
}

float particle::p(){
	return sqrtf(Px*Px + Py*Py + Pz*Pz);
}

float particle::aangle(){
	return asinf(Pz/p());
}

bool particle::of_type(int t, int c, bool m){
	return (mes == meson) && (c == charge) && (t == type);
}
