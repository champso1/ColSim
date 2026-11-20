#include "colsim/colsim.hpp"
using namespace colsim;

int main() {
	// create the main object
	ColSimMain colsim(ColSimMain::InitFlag::HARD_SCATTERING);

	// start/initialize generation
	// calculates cross section and other stuff
	// in preparation for event generation
	colsim.start();

	// generate events
	auto events = colsim.generate_events(1000000);

	// generate plots
	// colsim.generate_plots();

   	return 0;
}
