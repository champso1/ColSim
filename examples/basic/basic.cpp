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
	colsim.generateEvents(1000000);

	// generate plots
	colsim.generatePlots();

	// stop/deinitialize generation
	colsim.stop();
    
   	return 0;
}
