#include <iostream>

#include "particle.hpp"
#include "event.hpp"
#include "data.hpp"


int main(int argc, char **argv){

  if(argc < 4){
    puts("Usage: command <path to mesons file (fort.301)>"
	 " <path to baryons file (fort.300)>"
	 " <path to 'input' file>");
    return -1;
  }

  std::ifstream s(argv[3]);
  data D(s);
  s.close();
  // mesons
  {
    s.open(argv[1]);
    if(s.is_open()) D.readin_mesons(s);
    // else std::cout << "Warning: couldn't open mesons file "
    // 		   << argv[1] << '\n';
    s.close();
  }
  // baryons
  {
    s.open(argv[2]);
    if(s.is_open()) D.readin_baryons(s);
    // else std::cout << "Warning: couldn't open mesons file "
    // 		   << argv[2] << '\n';
    s.close();
  }

  // D.report_pnum(std::cout);

  std::cout << D.meson_count() << '\n';

  return 0;
}
