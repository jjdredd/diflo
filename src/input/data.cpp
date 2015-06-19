#include "data.hpp"


//
// particle
//

particle::particle(){}

particle::particle(double x, double y, double z, double Px,
		   double Py, double Pz, double P0, double b,
		   int type, int charge, bool meson)
	: x(x), y(y), z(z), Px(Px), Py(Py), Pz(Pz), P0(P0), b(b),
	  type(type), charge(charge), meson(meson){}

particle::~particle(){}

bool particle::operator==(double *p){
	return ((fabsf(this->P0 - p[0]) < THRES)
		&& (fabsf(this->Px - p[1]) < THRES)
		&& (fabsf(this->Py - p[2]) < THRES)
		&& (fabsf(this->Pz - p[3]) < THRES));
}

double particle::rapid(){
	return atanf(Py/P0);
}

double particle::p(){
	return sqrtf(Px*Px + Py*Py + Pz*Pz);
}

double particle::aangle(){
	return asinf(Pz/p());
}

bool particle::of_type(int t, int c, bool m){
	return (m == meson) && (c == charge) && (t == type);
}

//
// event
//

event::event(){}

event::~event(){}

void event::add_particle(particle &p){
	particles.push_back(p);
	return;
}

unsigned int event::particle_count(){
	return particles.size();
}

//
// data class
//

data::data(std::ifstream &s, HSDVersion v) : hsd_ver(v) {
	char str[128];
	int n;
	// read input
	s.getline(str, 128);
	sscanf(str, "%i", &A);
	for (n = 1; n < 5; n++) s.getline(str, 128); /* skip */
	sscanf(str, "%lf,", &Elab);
	for (; n < 9; n++) s.getline(str, 128);
	sscanf(str, "%i,", &NUM);
	s.getline(str, 128);
	sscanf(str, "%i", &ISUBS);
	s.getline(str, 128); 	/* skip seed */
	s.getline(str, 128);	/* time */
	sscanf(str, "%lf", &Time);

	P = new event* [ISUBS];
	for(int i = 0; i < ISUBS; i++) P[i] = new event[NUM];
}

data::~data(){
	for(int i = 0; i < ISUBS; i++)
		delete[] P[i];
	delete[] P;
}

bool data::parse_input_line(char *str, int *isub, int *irun, particle *p){

	switch(hsd_ver) {

	case HSD_VER_ORIG:
		return (sscanf(str, "\t%i\t%i\t%i\t%i\t%lf\t%lf\t%lf\t%lf"
			       "\t%lf\n", &p->type, &p->charge, isub, irun,
			       &p->Px, &p->Py, &p->Pz, &p->P0, &p->b) == 9);

	case HSD_VER_COORD:
		return (sscanf(str, "\t%i\t%i\t%i\t%i\t%lf\t%lf\t%lf\t%lf"
			       "\t%lf\t%lf\t%lf\t%lf\n",
			       &p->type, &p->charge, isub, irun,
			       &p->Px, &p->Py, &p->Pz, &p->P0,
			       &p->b, &p->x, &p->y, &p->z) == 12);

	case HSD_VER_PHSD:
	default:		// must be empty

		// TODO: for PHSD output format
		std::cerr << "HSD_VER_PHSD not implemented, exiting!\n";
		return false;
	}
}

void data::readin_particles(std::ifstream &s, bool mesons){
    
	int isub, irun;
	char str[256];
	particle ptcl;

	s.getline(str, 256);
	while(parse_input_line(str, &isub, &irun, &ptcl) && (!s.eof())){
		isub--; irun--;
		ptcl.meson = mesons;
		P[isub][irun].add_particle(ptcl);
		NParticles += 1.0;
		s.getline(str, 256);
	}
	return;
}

void data::report_pnum(std::ostream &os){
	int totnum = 0;
	for(int isub = 0; isub < ISUBS; isub++){
		for(int irun = 0; irun < NUM; irun++){
			totnum += P[isub][irun].particle_count();
		}
	}
	os << "Read in " << totnum << " particles\n";
	return;
}
