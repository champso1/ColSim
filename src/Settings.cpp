#include "ColSim/Settings.hpp"
#include "ColSim/Logger.hpp"
#include "ColSim/Utils.hpp"
#include "LHAPDF/Factories.h"

#include <cmath>
#include <string>
#include <vector>
#include <fstream>
#include <unordered_map>


namespace ColSim {

	Bool Settings::doesKeyExist(SettingsMap::const_iterator& it, const std::string& key) {
		Bool found;
		it = settings.find(key);
		found = (it != settings.end());

		return found;
	}
	

	void Settings::loadConfigFile(const std::string& configFilePath) {
		if (configFilePath.length() > 0) {
			std::ifstream configFileStream;
			configFileStream.open(configFilePath);
			if (!configFileStream) 
				LOGGER.logAbort("Could not open config file '%s'", configFilePath.c_str());

			std::string line;
			std::vector<std::string> tokens;
			while (std::getline(configFileStream, line)) {
				// skip empty lines and comments
				if ((line.size() <= 1) || (line[0] == '#'))
					continue;
		   
				SplitString(line, "=", tokens);

				settings.insert({tokens[0], tokens[1]});
			
				tokens.clear();
			}
			configFileStream.close();
		}

	    SettingsMap::const_iterator it;

		// center of mass energy
		// in TeV; must scale to GeV
		if (doesKeyExist(it, "ECM")) {
			ECM = std::stod(it->second);
			if (ECM < 5.0 || ECM > 20.0) {
				LOGGER.logAbort("ECM value of %.2lf is out of bounds of default [5.0,20.0] range."
								"Values outside this range are not recommended (LHC operates at 14TeV!).");
			}
			ECM *= 1000.0;
			S = ECM*ECM;
		} else {
			ECM = 14000.0;
			S = ECM*ECM;
		}
		LOGGER.logMessage("Using ECM=%lf", ECM);

		// name of PDF
		if(doesKeyExist(it, "PDFName"))
			PDFName = it->second;
		else
			PDFName = "cteq6l1";
		LHAPDF::setVerbosity(0);  // don't print anything, LHAPDF!
		pdf = LHAPDF::mkPDF(PDFName, 0);

		// process string
		Process = "PP2Zg2ll";
		LOGGER.logMessage("Process string=%s", Process.c_str());

		// number if iterations for cross section calculation: REQUIRED
		if(doesKeyExist(it, "NumXSIterations")) {
			numXSIterations = std::stoi(it->second);
			if (numXSIterations < 10000) {
				LOGGER.logWarning("%d Monte Carlo iterations is probably too low to ensure convergence."
								  "Perhaps try raising it.", numXSIterations);
			} else if (numXSIterations > 1e9) {
				LOGGER.logWarning("%e Monte carlo iterations is probably too high; the program will take a long time to run."
								  "Perhaps try lowering it.", numXSIterations);
			}
		} else
			numXSIterations = 1000000;
			
		
		LOGGER.logMessage("Using %d iterations for cross section calculation", numXSIterations);

		// cutoff energy for phase space generation
		if(doesKeyExist(it, "MinCutoffEnergy")) {
			minCutoffEnergy = std::stod(it->second);
			if (minCutoffEnergy < 1.0)
				LOGGER.logAbort("The specified cutoff energy %.2lf is less than the absolute minimum of 1.0.");
			else if (minCutoffEnergy > 500.0)
				LOGGER.logWarning("The specified cutoff energy %.2lf is higher than ordinary and will miss out on some phase space elements.");
			minCutoffEnergy_2 = minCutoffEnergy*minCutoffEnergy;
		} else {
			minCutoffEnergy = 60.0;
			minCutoffEnergy_2 = 60.0*60.0;
		}
		LOGGER.logMessage("Setting minimum cutoff energy for cross section calculation to %.3lf", minCutoffEnergy);

		// transformation mass/width (also for phase space points): REQUIRED
		if(doesKeyExist(it, "TransformationEnergy")) {
			transEnergy = std::stod(it->second);
			if (transEnergy < minCutoffEnergy) {
				LOGGER.logAbort("The specified transformation energy %.2lf cannot be less than Q_min.");
			} else if (minCutoffEnergy > std::sqrt(ECM))
				LOGGER.logAbort("The specified transformation energy %.2lf cannot be larger than sqrt(ECM).");
			transEnergy_2 = transEnergy*transEnergy;
		} else {
			transEnergy = minCutoffEnergy;
			transEnergy_2 = transEnergy*transEnergy;
		}
		LOGGER.logMessage("Setting transformation mass/energy to %.2lf", transEnergy);

		// initial evolution scale for parton showering: REQUIRED
		if(doesKeyExist(it, "InitialEvolEnergy")) {
			initialEvolEnergy = std::stod(it->second);
			if (initialEvolEnergy < 100.0 || initialEvolEnergy > 5000.0)
				LOGGER.logWarning("Initial evolution energy of %.2lf is out of ordinary range of [100.0,5000.0] GeV.");
			initialEvolEnergy_2 = initialEvolEnergy*initialEvolEnergy;
		} else {
			initialEvolEnergy = 1000.0;
			initialEvolEnergy_2 = 1000.0*1000.0;
		}
		LOGGER.logMessage("Will start parton evolution at %.2lf", initialEvolEnergy);


		// used fixed renormalization scale
		if(doesKeyExist(it, "FixedScale")) {
			fixedScale = settings.at("FixedScale").compare("Yes") == 0;
			if (fixedScale) {
				LOGGER.logMessage("Using fixed scale (mass of Z boson) for parton evolution.");
				LOGGER.logWarning("A fixed scale misses out on some higher order effects.");
			}
			else {
				LOGGER.logMessage("Using variable scale for parton evolution.");
			}
		}

		if(doesKeyExist(it, "EvolutionEnergyCutoff")) {
			evolEnergyCutoff = std::stod(it->second);
			if (evolEnergyCutoff < 1.0) 
				LOGGER.logAbort("Minimum cutoff energy for parton evolution MUST be greater than 1.0.");
			if (evolEnergyCutoff >= 100.0)
				LOGGER.logWarning("Specified cutoff energy for parton evolution is abnormally high.");
		} else {
			evolEnergyCutoff = 1.0;
		}
		LOGGER.logMessage("Will cut off parton evolution at %.3lf", evolEnergyCutoff);


		LOGGER.logMessage("-----------------------------------------");
	    LOGGER.logMessage("Finished loading configuration file data.");
		LOGGER.logMessage("-----------------------------------------");
	}




	void Settings::readString(const std::string& str) {
		std::vector<std::string> tokens;
		SplitString(str, "=", tokens);

		// if we have more than 2 tokens, print a warning and ignore this string
		if (tokens.size() != 2) {
			LOGGER.logWarning("(Settings::readString) invalid input string. perhaps there was an extra '=' sign?");
			return;
		}

		// simply iterate through all the keys
		const std::string& key = tokens.at(0);
		const std::string& val = tokens.at(1);

		if (key.compare("PDFName") == 0) {
			if (val.compare(PDFName) == 0) {
				LOGGER.logWarning("(Settings::readString) The PDF name '%s' is already in use.", val.c_str());
				return;
			}
			delete pdf;
			pdf = LHAPDF::mkPDF(val, 0);
			LOGGER.logMessage("Changing PDF from '%s' to '%s'", PDFName.c_str(), val.c_str());
			PDFName = val;
		} else if (key.compare("Process") == 0) {
			LOGGER.logWarning("(Settings::readString) there are no other possible processes. input string ignored");
			return;
		} else if (key.compare("NumXSIterations") == 0) {
			const UInt32 newNumXSIterations = std::stoi(val);
			if (newNumXSIterations == numXSIterations) {
				LOGGER.logWarning("(Settings::readString) The inputted number of Monte Carlo iterations is identical to what is already in use.");
				return;
			}
			LOGGER.logMessage("Changing number of Monte Carlo iterations from %d to %d", numXSIterations, newNumXSIterations);
			numXSIterations = newNumXSIterations;
		}else if (key.compare("ECM") == 0) {
			const Double newECM = std::stod(val)*1000.0, newS = newECM*newECM;
			if (newECM == ECM) {
				LOGGER.logWarning("(Settings::readString) The inputted ECM is identical to what is already in use.");
				return;
			}
			LOGGER.logMessage("Changing ECM from %.2lf to %.2lf (meaning S changes from %.2lf to %.2lf) TeV",
							  ECM/1000.0, newECM/1000.0, S/1000.0/1000.0, newS/1000.0/1000.0);
			ECM = newECM;
			S = newS;
		} else if (key.compare("MinCutoffEnergy") == 0) {
			const Double newMinCutoffEnergy = std::stod(val);
			if (newMinCutoffEnergy == minCutoffEnergy) {
				LOGGER.logWarning("(Settings::readString) The inputted minimum cutoff energy is identical to what is already in use.");
				return;
			}
			LOGGER.logMessage("Changing minimum cutoff energy from %.2lf to %.2lf",
							  minCutoffEnergy, newMinCutoffEnergy);
			minCutoffEnergy = newMinCutoffEnergy;
			minCutoffEnergy_2 = minCutoffEnergy*minCutoffEnergy;
		} else if (key.compare("TransformationEnergy") == 0) {
			const Double newTransEnergy = std::stod(val);
			if (newTransEnergy == transEnergy) {
				LOGGER.logWarning("(Settings::readString) The inputted transformation width/mass is identical to what is already in use.");
				return;
			}
			LOGGER.logMessage("Changing transformation mass/width from %.2lf to %.2lf",
							  transEnergy, newTransEnergy);
			transEnergy = newTransEnergy;
			transEnergy_2 = transEnergy*transEnergy;
		} else if (key.compare("InitialEvolEnergy") == 0) {
			const Double newInitialEvolEnergy = std::stod(val);
			if (newInitialEvolEnergy == initialEvolEnergy) {
				LOGGER.logWarning("(Settings::readString) The inputted initial evolution energy is identical to what is already in use.");
				return;
			}
			LOGGER.logMessage("Changing initial evolution energy from %.2lf to %.2lf",
							  initialEvolEnergy, newInitialEvolEnergy);
			initialEvolEnergy = newInitialEvolEnergy;
			initialEvolEnergy_2 = initialEvolEnergy*initialEvolEnergy;
		} else if (key.compare("FixedScale") == 0) {
			const Bool newFixedScale = val.compare("Yes") == 0;
			if (newFixedScale == fixedScale) {
				LOGGER.logWarning("(Settings::readString) The inputted fixed scale choice is identical to what is already in use.");
				return;
			}
			LOGGER.logMessage("Changing choice of fixed scale from %s to %s",
							  (fixedScale ? "true" : "false"),
							  (newFixedScale ? "true" : "false"));
			fixedScale = newFixedScale;
		} else if (key.compare("EvolutionEnergyCutoff") == 0) {
			const Double newEvolEnergyCutoff = std::stod(val);
			if (newEvolEnergyCutoff == evolEnergyCutoff) {
				LOGGER.logWarning("(Settings::readString) The inputted evolution energy cutoff is identical to what is already in use.");
				return;
			}
			LOGGER.logMessage("Changing evolution energy cutoff from %.2lf to %.2lf", evolEnergyCutoff, newEvolEnergyCutoff);
			evolEnergyCutoff = newEvolEnergyCutoff;
		} else {
			LOGGER.logWarning("(Settings::readString) unrecognized key. it it spelled correctly.");
		}
	}

}; // namespace ColSim
