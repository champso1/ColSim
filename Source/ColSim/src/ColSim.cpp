#include "ColSim/ColSim.hpp"

#include <fstream>
#include <stdexcept>

namespace ColSim {

	ColSimMain::ColSimMain(const std::string& configFilePath) {
		// initalize the logger with the default file path
	    LOGGER.initFile();

		// read in the settings from the configuration file
		SETTINGS.loadConfigFile(configFilePath);

		// load the hard process given in the config file
		loadHardProcess(SETTINGS.Process);
	}




	void ColSimMain::loadHardProcess(const std::string& processStr) {
		(void)processStr;
		hardProcess = std::unique_ptr<HardProcess>(new PP2Zg2ll());
	}


	void ColSimMain::CalcCrossSection(Double& result, Double& error) {
		MonteCarloParams params;
		params.num_evals = SETTINGS.numXSIterations;
		params.num_dims = hardProcess->phaseSpace->getDim();
		params.min = hardProcess->phaseSpace->getMins();
		params.max = hardProcess->phaseSpace->getMaxes();
		params.func = hardProcess->genDSigmaLambda();

		IntegrationResult res = MonteCarloIntegrate(params);
		result = res.result;
		error = res.error;
	}
}; // namespace ColSim
