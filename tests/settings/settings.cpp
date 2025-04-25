#include <ColSim/ColSim.hpp>
using namespace ColSim;

int main() {
	ColSimMain colsim;
	colsim.init(ColSimMain::HARD_SCATTERING);
	colsim.stop();	
   	return 0;
}
