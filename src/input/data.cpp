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
	return atanf(Pz/P0);
}

double particle::p(){
	return sqrtf(Px*Px + Py*Py + Pz*Pz);
}

double particle::aangle(){
	return asinf(Py/p());
}

bool particle::of_type(int t, int c, bool m){
	return (m == meson) && (c == charge) && (t == type);
}

int particle::octant(){
	if( p->Px > 0 ){
		if( p->Py > 0 ){
			if( p->Pz > 0 ) return 0;
			else return 1;
		}
		else{
			if( p->Pz > 0 ) return 2;
			else return 3;
		}
	}
	else{
		if( p->Py > 0 ){
			if( p->Pz > 0 ) return 4;
			else return 5;
		}
		else{
			if( p->Pz > 0 ) return 6;
			else return 7;
		}
	} 
}

double dotprod(particle &p1, particle &p2){
	return (p1.Px * p2.Px + p1.Py * p2.Py + p1.Pz * p2.Pz);
}

double dist(particle &p1, particle &p2){
	struct particle p;
	p.Px = p1.Px - p2.Px;
	p.Py = p1.Py - p2.Py;
	p.Pz = p1.Pz - p2.Pz;
	return p.p();
}

/* accepts sorted array */
double mixprod(particle **p){
	struct particle r;
	r.Px = p[1]->Py*p[0]->Pz - p[1]->Pz*p[0]->Py;
	r.Py = p[1]->Pz*p[0]->Px - p[1]->Px*p[0]->Pz;
	r.Pz = p[1]->Px*p[0]->Py - p[1]->Py*p[0]->Px;
	return dotprod(p[2], &r);
}

int pcompare(particle &a, particle &b){
	double c = dotprod(a, a) - dotprod(b, b);
	if(c > 0) return 1;
	else{
		if(c < 0) return -1;
		else return 0;
	}
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

data::data(std::ifstream &s, DataVersion v) : hsd_ver(v) {
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
	for(unsigned i = 0; i < ISUBS; i++) P[i] = new event[NUM];
}

data::~data(){
	for(unsigned i = 0; i < ISUBS; i++)
		delete[] P[i];
	delete[] P;
}

bool data::parse_input_line(char *str, int *isub, int *irun, particle *p){

	int unk;

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
		return (sscanf(str, "\t%i\t%i\t%i\t%i\t%lf\t%lf\t%lf\t%lf"
			       "\t%lf\t%i\n",
			       &p->type, &p->charge, isub, irun,
			       &p->Px, &p->Py, &p->Pz, &p->P0,
			       &p->b, &unk) == 10 );

	default:		// must be empty

		// TODO: for PHSD output format
		std::cerr << "DataVersion not implemented, exiting!\n";
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
		NParticles += 1;
		s.getline(str, 256);
	}
	return;
}

unsigned data::NumberOfParticles(){
	unsigned totnum = 0;
	for(unsigned isub = 0; isub < ISUBS; isub++){
		for(unsigned irun = 0; irun < NUM; irun++){
			totnum += P[isub][irun].particle_count();
		}
	}
	return totnum;
}


//
// DataSIn class
//


DataSIn::DataSIn() : EventNum(0), isub(0), irun(0) {}

DataSIn::~DataSIn(){}

int DataSIn::parse_input_line(char *str, particle *p){

	int unk;

	switch(dat_ver) {

	case HSD_VER_ORIG:
		if (sscanf(str, "\t%i\t%i\t%i\t%i\t%lf\t%lf\t%lf\t%lf"
			   "\t%lf\n", &p->type, &p->charge, isub, irun,
			   &p->Px, &p->Py, &p->Pz, &p->P0, &p->b) == 9){}

	case HSD_VER_COORD:
		if (sscanf(str, "\t%i\t%i\t%i\t%i\t%lf\t%lf\t%lf\t%lf"
			   "\t%lf\t%lf\t%lf\t%lf\n",
			   &p->type, &p->charge, isub, irun,
			   &p->Px, &p->Py, &p->Pz, &p->P0,
			   &p->b, &p->x, &p->y, &p->z) == 12){

			if (isub != sub || irun != run) return EventNum + 1;
			else return EventNum;

		} else return -1;

	case HSD_VER_PHSD:
		if (sscanf(str, "\t%i\t%i\t%i\t%i\t%lf\t%lf\t%lf\t%lf"
			   "\t%lf\t%i\n",
			   &p->type, &p->charge, isub, irun,
			   &p->Px, &p->Py, &p->Pz, &p->P0,
			   &p->b, &unk) == 10 ){}

	case ALICE_GUYS:

	case ROGACH:

	default:		// must be empty

		// TODO: for PHSD output format
		std::cerr << "DataVersion not implemented, exiting!\n";
		return false;
	}

}

bool DataSIn::FetchEvent(std::ifstream &s, event &e){

}

void DataSIn::readin_data(std::ifstream &s){

}
