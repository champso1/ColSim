#include "ColSim/ColSim.hpp"
using namespace ColSim;
using namespace std;

#include <iostream>

int main() {
	// create the main object
	ColSimMain colsim;

	// start/initialize generation
	colsim.start();

	// generate 10 events
	for (int i=0; i<100; i++)
		cout << colsim.generateEvent() << "\n";

	// stop/deinitialize generation
	colsim.stop();
    
   	return 0;
}
