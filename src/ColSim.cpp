#include "ColSim/ColSim.hpp"

#include "ColSim/Constants.hpp"
#include "ColSim/Logger.hpp"
#include "ColSim/PartonShower.hpp"
#include "ColSim/Settings.hpp"
#include "ColSim/PhaseSpace.hpp"
#include "ColSim/Math.hpp"
#include "ColSim/Utils.hpp"
#include "ColSim/Gnuplot.hpp"

#include <algorithm>
#include <fstream>
#include <memory>
#include <chrono>
#include <stdexcept>

namespace ColSim {

	ColSimMain::ColSimMain(const std::string& logFilePath) {
	    // initalize the logger with the default file path
	    LOGGER.initFile(logFilePath);

 		// -- random number generation --
		UInt32 seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
		engine.seed(seed);
		//engine.seed(DEFAULT_RANDOM_SEED);
		randDouble = std::bind(distribution, engine);
	}

	ColSimMain::~ColSimMain() {
		LOGGER.logMessage("Shutting down...");
	}


	void ColSimMain::init(InitFlag initFlag, const std::string& configFilePath) {
		flag = initFlag;
		
		// read in the settings from the given configuration file
		SETTINGS.loadConfigFile(configFilePath);

		// load the process given in the config file
		switch(initFlag) {
			case HARD_SCATTERING:
				loadHardProcess(SETTINGS.Process);
				break;
			case PARTON_SHOWERING:
				partonShower = std::unique_ptr<PartonShower>(new GluonShower());
				break;
		}
	}




	void ColSimMain::loadHardProcess(const std::string& processStr) {
		if (processStr.compare("PP2Zg2ll") != 0)
			LOGGER.logAbort("Only 'PP2Zg2ll' has been implemented.");

		hardProcess = std::unique_ptr<HardProcess>(new PP2Zg2ll());
		//hardProcess = std::unique_ptr<HardProcess>(new PP2Jets());
	}


	void ColSimMain::start() {
		switch(flag) {
			case HARD_SCATTERING:
				start_hardProcess();
				break;
			case PARTON_SHOWERING:
				start_partonShower();
				break;
		}
	};
	

	Bool ColSimMain::generateEvent() {
		switch(flag) {
			case HARD_SCATTERING:
				return generateEvent_hardProcess();
			case PARTON_SHOWERING:
				return generateEvent_partonShower();
		}

		return true; // unreachable, but gcc shuts the fuck up
	}

	void ColSimMain::generateEvents(UInt32 numEvents) {
		plotPoints.clear();
		emissionRecord.clear();
		while (numEvents > 0) {
			if(generateEvent())
				numEvents--;
		}
	}


	void ColSimMain::generatePlots() {
		switch(flag) {
			case HARD_SCATTERING:
				generatePlots_hardProcess();
				break;
			case PARTON_SHOWERING:
				generatePlots_partonShower();
				break;
		}
	}

	
	void ColSimMain::stop() {
		// does nothing at the moment!
		LOGGER.logMessage("Stopped event generation!");
		return;
	}




	Bool ColSimMain::generateEvent_hardProcess() {
		PhaseSpace phaseSpace = hardProcess->getPhaseSpace();

		std::vector<Double> phaseSpacePoints;
		phaseSpace.fillPhaseSpace(phaseSpacePoints);
		const std::vector<Double>& deltas = phaseSpace.getDeltas();
		
		HardProcess::Result res = hardProcess->dSigma(phaseSpacePoints);
		if (res == HardProcess::Result::invalidResult())
			return false;
		
		Double weight = res.weight;
		// scale the weight
		for (Double d : deltas)
			weight *= d;

		// hit-or-miss to see if it actually works
		Double weightRatio = weight/maxWeight;
		Double randNum = randDouble();

		while (randNum > weightRatio) {
			phaseSpace.fillPhaseSpace(phaseSpacePoints);
			res = hardProcess->dSigma(phaseSpacePoints);
			weight = res.weight;
			for (Double d : deltas)
				weight *= d;

			// ensure to calculate the new weight and another random number
			weightRatio = weight/maxWeight;
			randNum = randDouble();
		}

		// generate a list of particles
		std::vector<Particle> particles;
		hardProcess->generateParticles(particles);
		eventRecord.emplace_back(Event(weight, particles));

		// plot the phase space points and the chosen additional values
		std::vector<Double> _plotPoints = Join(phaseSpacePoints, res.additionalVals);
	    plotPoints.emplace_back(_plotPoints);

		return true;
	}

	Bool ColSimMain::generateEvent_partonShower() {
	    const Double Q0 = SETTINGS.initialEvolEnergy;
		const Double Qf = SETTINGS.evolEnergyCutoff;
		
		Double scale;
		if (SETTINGS.fixedScale)
			scale = Q0/2.0;
		else {
			scale = Qf;
		}

		const AlphaS& asRef = partonShower->getAlphaSRef();
		emissionRecord.emplace_back(partonShower->Evolve(Q0, Qf, asRef.getAlphaSOver(scale)));
		return true;
	}


	void ColSimMain::start_hardProcess() {
		// ensure to reset ranges in the event that
		// we have changed a config file variable between calculations
		hardProcess->getPhaseSpace().setRanges();
		
		LOGGER.logMessage("Calculating cross section via Monte Carlo integration...");
		HardProcessResult res = hardProcess->calculate();
		LOGGER.logMessage("Finished!");
		LOGGER.logMessage("Result is: %.9lf +- %.9lf pb (picobarns)",
						  res.result, res.error);

		crossSection = res.result;
		crossSectionError = res.error;
		maxWeight = res.maxWeight;
		maxPSPoints = res.maxPoints;

		LOGGER.logMessage("Maximum weight achieved: %.9lf", maxWeight);
	}

	void ColSimMain::start_partonShower() {
		// nothing else to initialize
		LOGGER.logMessage("Starting parton shower event generation.");
	}


	void ColSimMain::generatePlots_hardProcess() {
		LOGGER.logMessage("Generating plots...");

		PhaseSpace& phaseSpace = hardProcess->getPhaseSpace();

		Gnuplot plot;
		plot.setHistInfo(phaseSpace.getMins(),
						 phaseSpace.getMaxes(),
						 phaseSpace.getDeltas(), 100);
		std::vector<std::string> plotColNames(phaseSpace.getNames());
		plot.openDataFile("events.dat", plotColNames);
	    plot.addDataPoints(plotPoints);

		LOGGER.logMessage("Saved datafile in 'events.dat'");
		LOGGER.logMessage("Creating temporary Gnuplot script files for plot generation...");

		// plot.setTitles(phaseSpace.getTitles());
		plot.setXLabels(phaseSpace.getXLabels());
		plot.setYLabels(phaseSpace.getYLabels());
		
		plot.plot();

		LOGGER.logMessage("Plots saved! Check the 'plots' directory to view them.");
	}

	
	void ColSimMain::generatePlots_partonShower() {
		const std::vector<std::string> plotNames{"t", "pT", "m"};

		std::vector<Double> v;
	    for (const std::vector<PartonShower::Emission>& emissions : emissionRecord) {
			for (const PartonShower::Emission& e : emissions) {
				v.push_back(std::sqrt(e.t));
				v.push_back(std::sqrt(e.pT_2));
				v.push_back(std::sqrt(e.m_2));
				plotPoints.emplace_back(v);
				v.clear();
			}
		}

		// maximum for t is directly cutoff at this value,
		// so we can fix it here
		double tMax = std::sqrt(SETTINGS.initialEvolEnergy);
		
		std::vector<Double>
			min{0.0, 0.0, 0.0},
			max{tMax, 20.0, 25.0},
			delta{tMax, 20.0, 25.0};

		
		Gnuplot plot;
		plot.setHistInfo(min, max, delta, 100);
		plot.setXLabels({"√t", "p_T", "m_{virt}"});
		plot.setYLabels({"Events", "Events", "Events"});
		plot.openDataFile("emissions.dat", plotNames);
		plot.addDataPoints(plotPoints);

		LOGGER.logMessage("Saved datafile in 'emissions.dat'");
		LOGGER.logMessage("Creating temporary Gnuplot script files for plot generation...");

		
		
		plot.plot();

		LOGGER.logMessage("Plots saved! Check the 'plots' directory to view them.");
	}
}; // namespace ColSim
