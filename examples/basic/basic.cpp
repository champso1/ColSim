#include "ColSim/ColSim.hpp"
using namespace ColSim;
using namespace std;

#include <iostream>
#include <fstream>

int main() {
	// create the main object
	ColSimMain colsim;

	colsim.init(ColSim::ColSimMain::HARD_SCATTERING);

	// start/initialize generation
	colsim.start();

	// generate events
	colsim.generateEvents(100000);

	// stop/deinitialize generation
	colsim.stop();

	std::vector<double> pt, e;
	for (const Event& event: colsim.getEventRecord()) {
		
	}

	std::ofstream outfile("data.dat");
	for (uint i=0; i<pt.size(); i++) {
		outfile << pt.at(i) << '\t' << e.at(i) << '\n';
	}
	outfile.close();
    
   	return 0;
}
