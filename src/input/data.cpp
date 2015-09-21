#include "data.hpp"


//
// particle
//

particle::particle() {}

particle::particle(double x, double y, double z, double Px,
		   double Py, double Pz, double P0, double b,
		   int type, int charge, bool meson)
	: x(x), y(y), z(z), Px(Px), Py(Py), Pz(Pz), P0(P0), b(b),
	  type(type), charge(charge), meson(meson) {}

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

unsigned particle::octant(){
	if( Px > 0 ){
		if( Py > 0 ){
			if( Pz > 0 ) return 0;
			else return 1;
		}
		else{
			if( Pz > 0 ) return 2;
			else return 3;
		}
	}
	else{
		if( Py > 0 ){
			if( Pz > 0 ) return 4;
			else return 5;
		}
		else{
			if( Pz > 0 ) return 6;
			else return 7;
		}
	} 
}

unsigned particle::diant() {
	if ((-Px*tan(RPAngle) + Py) > 0) return 0;
	else return 1;
}

double dotprod(const particle &p1, const particle &p2){
	return (p1.Px * p2.Px + p1.Py * p2.Py + p1.Pz * p2.Pz);
}

double dist(const particle &p1, const particle &p2){
	struct particle p;
	p.Px = p1.Px - p2.Px;
	p.Py = p1.Py - p2.Py;
	p.Pz = p1.Pz - p2.Pz;
	return p.p();
}

/* accepts sorted array */
double mixprod(const particle &a, const particle &b, const particle &c ){
	struct particle r;
	r.Px = b.Py * a.Pz - b.Pz * a.Py;
	r.Py = b.Pz * a.Px - b.Px * a.Pz;
	r.Pz = b.Px * a.Py - b.Py * a.Px;
	return dotprod(c, r);
}

bool pcompare(const particle a, const particle b){
	return dotprod(a, a) < dotprod(b, b);
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

data::data(std::ifstream &s, DataVersion v, bool pick = false, int type = 0)
	: hsd_ver(v), pick(pick), type(type) {
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
		if (pick && type == ptcl.type) {
			ptcl.meson = mesons;
			P[isub][irun].add_particle(ptcl);
			NParticles += 1;
		}
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
// ALICE event fetcher
//

ALICEData::ALICEData(const char *file) : s(file) {}

ALICEData::~ALICEData(){}

bool ALICEData::parse_hline(char *str) {
	return (sscanf(str, "%i\t%i\n", &cur_nev, &cur_npart) == 2);
}

bool ALICEData::parse_line(char *str, particle &p) {
	int i;
	return (sscanf(str, "%i\t%lf\t%lf\t%lf\n", &i, &p.Px, &p.Py, &p.Pz)
		== 4);
}

bool ALICEData::FetchEvent(event &e){
	particle p;
	char str[256];
	s.getline(str, 256);
	if (!parse_hline(str)) return false;
	e.particles.clear();
	e.particles.reserve(cur_npart);
	for (unsigned i = 0; i < cur_npart; i++) {
		s.getline(str, 256);
		if (!parse_line(str, p)) return false;
		e.add_particle(p);
	}
	return true;
}

// fetch an event that has at least 'lower' particles and at most 'upper'
// particles. If upper is '0', the upper limit is removed
bool ALICEData::FetchNumEvent(event &e, unsigned lower, unsigned upper = 0) {
	if (upper && lower > upper) return false;
	while (FetchEvent(e)) {
		if (cur_npart > lower
		    && (!upper || cur_npart < upper))
			return true;
	}
	return false;
}


#if 0

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

#endif	// comment out DataSIn class
