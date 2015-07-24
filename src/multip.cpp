#include <iostream>
#include <fstream>

#include "data.hpp"

bool is_strange(particle &p) {

}

int main(int argc, char **argv){

	if(argc < 4){
		puts("Usage: command <path to mesons file (fort.301)>"
		     " <path to baryons file (fort.300)>"
		     " <path to 'input' file>");
		return -1;
	}

	//
	// TODO getopt for advanced tweaking
	//

	std::ifstream s(argv[3]);
	data D(s, HSD_VER_PHSD);
	s.close();
	// mesons
	{
		s.open(argv[1]);
		if(s.is_open()) D.readin_particles(s, true);
		// else std::cout << "Warning: couldn't open mesons file "
		// 		   << argv[1] << '\n';
		s.close();
	}
	// baryons
	{
		s.open(argv[2]);
		if(s.is_open()) D.readin_particles(s, false);
		// else std::cout << "Warning: couldn't open mesons file "
		// 		   << argv[2] << '\n';
		s.close();
	}

	std::cout << "Read in " << D.NumberOfParticles() << " particles"
		  << std::endl;

	

	return 0;
}
