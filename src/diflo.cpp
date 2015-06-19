#include <iostream>
#include <fstream>

#include "dist.hpp"
#include "data.hpp"
#include "momentum.hpp"


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
	data D(s, HSD_VER_COORD);
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

	D.report_pnum(std::cout);

	std::cout << "<P_x>, for P_z > 0 : " << MeanPx(D, true)
		  << "\n<P_x>, for P_z < 0 : " << MeanPx(D, false) << "\n";

	distribution dst(25, 100);
	dst.DataDist(D);
	dst.DistTransform(0.5);
	dst.PrintFlows(std::cout);

	std::ofstream distout("./distout");
	dst.WriteDistr(distout);
	std::ofstream ftout("./ftout");
	dst.WriteDistrFT(ftout);

	return 0;
}
