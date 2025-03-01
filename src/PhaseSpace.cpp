#include "ColSim/PhaseSpace.hpp"

#include "ColSim//Settings.hpp"

#include <random>

namespace ColSim {
	PhaseSpace::~PhaseSpace() {
		delete[] min;
		delete[] max;
		delete[] delta;
	}
	
	PhaseSpace_TauYCosth::PhaseSpace_TauYCosth() : PhaseSpace(3) {
		// initialize rho ranges
		const Double S = SETTINGS.S;
		rhoMin = std::atan((Q_MIN_2-MASS_TR_2) / (WIDTH_TR*MASS_TR));
		rhoMax = std::atan((S-MASS_TR_2) / (WIDTH_TR*MASS_TR));
		deltaRho = rhoMax - rhoMin;

		// set ranges
		min[0] = -1.0; max[0] = 1.0;      // cosTheta
		min[1] = rhoMin; max[1] = rhoMax; // rho
		min[2] = 0.0; max[2] = 1.0;       // y
	}





	PhaseSpace_Costh::PhaseSpace_Costh() : PhaseSpace(3) {
		*min = -1.0; *max = 1.0;
	}
}; // namespace ColSim
