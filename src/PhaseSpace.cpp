#include "ColSim/PhaseSpace.hpp"

#include "ColSim//Settings.hpp"
#include "ColSim/Math.hpp"

#include <random>

namespace ColSim {
	void PhaseSpace::fillDelta() {
		delta.clear();
		for (UInt i=0; i<min.size(); i++)
			delta.push_back(max.at(i) - min.at(i));
	}

	void PhaseSpace::fillPhaseSpace(std::vector<Double>& vec) {
		vec.clear();
	    for (UInt32 i=0; i<numDims; i++) {
			Double val = randDouble()*delta[i] + min[i];
			vec.push_back(val);
		}
	}

	
	
	
	PhaseSpace_TauYCosth::PhaseSpace_TauYCosth()
		: PhaseSpace(3,
					 {"cos(theta)", "rho", "y", "Q", "x1"},
					 {"cos(θ)", "ρ", "rapidity", "Q", "x_1"},
					 {"cos(θ)", "ρ", "y", "Q [GeV]", "x_1"},
					 {"Events", "Events", "Events", "Events", "Events"})
	{
		
	    setRanges();
	}

	void PhaseSpace_TauYCosth::setRanges() {
		const Double S = SETTINGS.S;
		const double E_TR_2 = SETTINGS.transEnergy_2;
		const double Q_MIN_2 = SETTINGS.minCutoffEnergy_2;

		// initialize rho ranges
		// rhoMin = std::atan((Q_MIN_2-E_TR_2) / (E_TR_2));
		// rhoMax = std::atan((S-E_TR_2) / (E_TR_2));
		rhoMin = std::atan(Q_MIN_2-E_TR_2) / E_TR_2;
		rhoMax = std::atan((S-E_TR_2) / (E_TR_2));

		
		
		deltaRho = rhoMax - rhoMin;

		min.clear(); max.clear();
		// set ranges for main variables
		min.push_back(-1.0); max.push_back(1.0);      // cosTheta
		min.push_back(rhoMin); max.push_back(rhoMax); // rho
		min.push_back(0.0); max.push_back(1.0);       // y

		// Q and x1, which aren't themselves
		// independent phase space variables
		// but are variables of interest later on
		min.push_back(60.0); max.push_back(500.0);  // Q
		min.push_back(0.0); max.push_back(1.0);    // x

		// fill the deltas
		fillDelta();
	}







	PhaseSpace_EtEta::PhaseSpace_EtEta()
		: PhaseSpace(2, {"E_t", "eta"})
	{
		setRanges();
	}


	void PhaseSpace_EtEta::setRanges() {
		// set ranges for main variables
		min.push_back(0.1); max.push_back(500.0); // E_t
		min.push_back(-5.0); max.push_back(5.0);    // eta

		// fill the deltas
		fillDelta();
	}
}; // namespace ColSim
