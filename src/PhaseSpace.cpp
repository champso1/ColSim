#include "ColSim/PhaseSpace.hpp"

#include "ColSim//Settings.hpp"
#include "ColSim/Math.hpp"

#include <random>

namespace ColSim {
	void PhaseSpace::fillDelta() {
		for (int i=0; i<getDim(); i++)
			delta[i] = max[i] - min[i];
	}

	void PhaseSpace::fillPhaseSpace(std::vector<Double>& vec) {
		// ensure that there is capacity to do this
		vec.resize(3);
	    for (UInt32 i=0; i<getDim(); i++) {
			Double val = randDouble()*getDeltas()[i] + getMins()[i];
			vec[i] = val;
		}
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
		fillDelta();
	}
}; // namespace ColSim
