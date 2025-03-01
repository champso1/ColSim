#include "ColSim/ColSim.hpp"
using namespace ColSim;

#include <iostream>

int main() {

	ColSimMain colsim;
	Double result, error;
	colsim.CalcCrossSection(result, error);
	std::cout << result << " +- " << error << std::endl;
	return 0;
}
