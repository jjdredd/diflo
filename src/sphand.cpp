#include <iostream>
#include <string>
#include <boost/program_options.hpp>

#include "data.hpp"
#include "grid.hpp"

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

	std::cout << cmdvars["dir"].as<std::string>()
		  << '\t' << cmdvars["out"].as<std::string>() << std::endl;

	std::ifstream s(cmdvars["dir"].as<std::string>() + std::string("input"));
	DataHSD D(s, false, 2);
	s.close();

	
	

	return 0;
}
