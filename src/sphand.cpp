#include <iostream>
#include <boost/program_options.hpp>

#include "input.hpp"
#include "grid.hpp"

namespace po = boost::program_options;

po::variables_map parse_options(int argc, const char **argv) {
	po::variables_map clvars;
	po::options_description desc("options");

	desc.add_options()
		("help,h", "Short Help description")
		("dir,d", "directory with Input file for hsd")
		("out,o", "output basename");

	try{
		po::store(po::parse_command_line(argc, argv, desc), clvars);
	} catch (po::error &e) {
		std::cerr << e.what() << std::endl;
	}
	return clvars;
}

int main(int argc, const char **argv) {
	po::variables_map
}
