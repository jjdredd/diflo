#include <iostream>
#include <string>
#include <boost/program_options.hpp>

#include "data.hpp"
#include "grid.hpp"
#include "handedness.hpp"
#include "helicity.hpp"

namespace po = boost::program_options;

po::variables_map parse_options(int argc, const char **argv) {
	po::variables_map clvars;
	po::options_description desc("Options");

	desc.add_options()
		("help,h", "This help description")

		("dir,d", po::value<std::string>()->required(),
		 "directory with Input file for hsd")

		("out,o", po::value<std::string>()->required(),
		 "output basename");

	try {
		po::store(po::parse_command_line(argc, argv, desc), clvars);
	} catch (po::required_option &e) {
		std::cerr << "Missing required option" << std::endl;
		std::cerr << e.what() << std::endl;
		exit(-1);
	} catch (po::error &e) {
		std::cerr << e.what() << std::endl;
		exit(-1);
	}

	if (clvars.count("help")) {
		std::cout << desc << std::endl;
	}
	po::notify(clvars);

	return clvars;
}

int main(int argc, const char **argv) {
	po::variables_map cmdvars = parse_options(argc, argv);

	if (cmdvars.count("help")) {
		return 0;
	}

	std::ifstream s(cmdvars["dir"].as<std::string>() + std::string("input"));
	DataHSDC D(s, false, 2);
	s.close();

	{
		s.open(cmdvars["dir"].as<std::string>() + std::string("fort.301"));
		if(s.is_open()) D.readin_particles(s, true);
		else std::cout << "Warning: couldn't open mesons file " << std::endl;
		s.close();
	}

	{
		s.open(cmdvars["dir"].as<std::string>() + std::string("fort.300"));
		if(s.is_open()) D.readin_particles(s, true);
		else std::cout << "Warning: couldn't open mesons file " << std::endl;
		s.close();
	}

	SymGrid g(12, 50);
	ParticleGrid pg(g, 4);

	for(unsigned isub = 0; isub < D.ISUBS; isub++){
		for(unsigned irun = 0; irun < D.NUM; irun++){
			pg.Populate(D.P[isub][irun]);
		}
	}

	// for(unsigned isub = 0; isub < D.ISUBS; isub++){
	// 	pg.Populate(D.P[isub][0]);
	// }

	pg.ShrinkToFit();

	pg.WriteParticleCount(std::string(cmdvars["out"].as<std::string>())
			      + std::string("pcnt_"));

	HandednessGrid HG(g);
	HG.Compute(pg);
	HG.WriteOutHandedness(std::string(cmdvars["out"].as<std::string>())
			      + std::string("handgrid_"));

	Helicity Hel(g);
	Hel.Compute(pg);
	Hel.WriteOutHelicity(std::string(cmdvars["out"].as<std::string>())
			      + std::string("helgrid_"));

	return 0;
}
