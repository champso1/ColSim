#include "ColSim/ColSim.hpp"

#include "ColSim/Constants.hpp"
#include "ColSim/Logger.hpp"
#include "ColSim/PartonShower.hpp"
#include "ColSim/Settings.hpp"

#include "ColSim/PhaseSpace.hpp"
#include "ColSim/Math.hpp"
#include "ColSim/Utils.hpp"

#include <algorithm>
#include <fstream>
#include <memory>
#include <chrono>
#include <stdexcept>

namespace ColSim {

	ColSimMain::ColSimMain(const std::string& configFilePath, const std::string& logFilePath) {
	    // initalize the logger with the default file path
	    LOGGER.initFile(logFilePath);

		// read in the settings from the configuration file
		SETTINGS.loadConfigFile(configFilePath);

		// load the hard process given in the config file
		loadHardProcess(SETTINGS.Process);
		loadPartonShower(SETTINGS.doPhotonEmission, SETTINGS.doGluonEmission);


		// other miscellaneous setup

		// -- random number generation --
		//UInt32 seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
		//engine.seed(seed);
		engine.seed(DEFAULT_RANDOM_SEED);
		randDouble = std::bind(distribution, engine);
	}




	void ColSimMain::loadHardProcess(const std::string& processStr) {
		std::vector<std::string> processTokens;
		SplitString(processStr, "2", processTokens);

		// ensure there are three "steps", i.e. an initial state,
		// an intermediate gauge boson, and a final state
		if (processTokens.size() == 2) {
			LOGGER.logAbort("Must specify an intermediate particle!");
		} else if ((processTokens.size() < 2)
				   || (processTokens.size() > 3)) {
			LOGGER.logAbort("Invalid number of process steps.");
		}

		// determine if the initial state is PP or ee
		// (anything else is forbidden currently)
		if ((processTokens[0].compare("PP") != 0)
			&& (processTokens[0].compare("ee") != 0)) {
			LOGGER.logAbort("Initial state must be PP or ee.");
		}

		// determine if the intermediate state is Zg or not
		if (processTokens[1].compare("Zg") != 0) {
			LOGGER.logAbort("Intermediate state must be Zg");
		}

		// lastly, ensure that the final state is two leptons
		std::vector<std::string> availableFinalStates{
			"ee", "ll", "mumu"
		};
		if (std::find(availableFinalStates.begin(), availableFinalStates.end(),
					  processTokens[2]) == availableFinalStates.end()) {
			LOGGER.logAbort("Invalid final state. Must be ee, ll or mumu.");
		}

		// now set the hard process according to the final state
		if (processTokens[0].compare("ee") == 0)
			LOGGER.logAbort("ee initial state not implemented!");
		else {
			if (processTokens[2].compare("ee") == 0)
				LOGGER.logAbort("ee final state not implemented!");
			else if (processTokens[2].compare("ll") == 0) {
				hardProcess = std::unique_ptr<HardProcess>(new PP2Zg2ll());
				return;
			} else
				LOGGER.logAbort("mumu final state not implemented");
		}
		
		LOGGER.logAbort("UNREACHABLE");
	}

	void ColSimMain::loadPartonShower(Bool doPhotonEmission, Bool doGluonEmission) {
		if (doPhotonEmission) {}
		else if (doGluonEmission)
			partonShower = std::unique_ptr<PartonShower>(new GluonShower());
		else
		    LOGGER.logError("Invalid combination of photon/gluon emission bits.");
	}

	void ColSimMain::start() {
		LOGGER.logMessage("Calculating cross section via Monte Carlo integration...");
		HardProcessResult res = hardProcess->calculate();
		LOGGER.logMessage("Finished!");
		LOGGER.logMessage("Result is: %.9lf +- %.9lf pb (picobarns)",
						  res.result, res.error);

		crossSection = res.result;
		crossSectionError = res.error;
		maxWeight = res.maxWeight;
		maxPSPoints = std::move(res.maxPoints); // not that it matters

		// log some of this stuff
		LOGGER.logMessage("Maximum weight achieved: %.9lf", maxWeight);
	};
	

	const Event& ColSimMain::generateEvent() {
		// todo: make this actually do something!
		
		// std::vector<PartonShower::Emission> emissions = partonShower->Evolve(1000.0, 1.0, 0.01627095);

		std::vector<Double> phaseSpacePoints, deltas;
		hardProcess->getPhaseSpace().fillPhaseSpace(phaseSpacePoints);
		deltas = hardProcess->getPhaseSpace().getDeltas();
		
		Double weight = hardProcess->dSigma(phaseSpacePoints);
		for (Double d : deltas)
			weight *= d;

		// hit-or-miss to see if it actually works
		Double weightRatio = weight/maxWeight;
		Double randNum = randDouble();

		UInt32 numRejectedEvents = 0;
		
		while (randNum > weightRatio) {
			numRejectedEvents++;
			hardProcess->getPhaseSpace().fillPhaseSpace(phaseSpacePoints);
			weight = hardProcess->dSigma(phaseSpacePoints);
			for (Double d : deltas)
				weight *= d;

			weightRatio = weight/maxWeight;
			randNum = randDouble();
		}
		LOGGER.logMessage("Rejected %u events for this iteration.",
						  numRejectedEvents);
		
		std::vector<Particle> particles;
		hardProcess->generateParticles(particles);

		eventRecord.emplace_back(Event(weight, particles));
		
		return eventRecord.back();
	}

	
	void ColSimMain::stop() {
		// does nothing at the moment!
		return;
	}
}; // namespace ColSim
