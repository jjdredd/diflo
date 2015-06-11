#include "data.hpp"

data::data(std::ifstream &s){
	char str[128];
	int n;
	// read input
	s.getline(str, 128);
	sscanf(str, "%i", &A);
	for (n = 1; n < 5; n++) s.getline(str, 128); /* skip */
	sscanf(str, "%f,", &Elab);
	for (; n < 9; n++) s.getline(str, 128);
	sscanf(str, "%i,", &NUM);
	s.getline(str, 128);
	sscanf(str, "%i", &ISUBS);
	s.getline(str, 128); 	/* skip seed */
	s.getline(str, 128);	/* time */
	sscanf(str, "%f", &Time);

	P = new event* [ISUBS];
	for(int i = 0; i < ISUBS; i++) P[i] = new event[NUM];
}

data::~data(){
	for(int i = 0; i < ISUBS; i++)
		delete[] P[i];
	delete[] P;
}

void data::readin_particles(std::ifstream &s, bool mesons){
    
	int isub, irun;
	char str[256];
	particle ptcl;

	s.getline(str, 256);
	while((sscanf(str, "\t%i\t%i\t%i\t%i\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",
		      &ptcl.type, &ptcl.charge, &isub, &irun,
		      &ptcl.Px, &ptcl.Py, &ptcl.Pz, &ptcl.P0, &ptcl.b,
		      &ptcl.x, &ptcl.y, &ptcl.z) == 12)
	      && (!s.eof())){

		isub--; irun--;
		ptcl.meson = mesons;
		P[isub][irun].add_particle(ptcl);
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
