#include "colsim/colsim.hpp"
using namespace colsim;
using namespace std;

int main() {
	// create the main object
	ColSimMain colsim;

	colsim.init(ColSimMain::HARD_SCATTERING);

	// start/initialize generation
	colsim.start();

	// generate events
	colsim.generate_events(1000000);

	// generate plots
	// colsim.generatePlots();

	// stop/deinitialize generation
	colsim.stop();
    
   	return 0;
}
