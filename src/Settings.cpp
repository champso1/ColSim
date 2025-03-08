#include "ColSim/Settings.hpp"
#include "ColSim/Logger.hpp"
#include "ColSim/Utils.hpp"

#include <cmath>
#include <vector>
#include <fstream>
#include <unordered_map>


namespace ColSim {

	void Settings::loadConfigFile(const std::string& configFilePath) {
		typedef std::unordered_map<std::string, std::string> SettingsMap;
		
		std::ifstream configFileStream;
		configFileStream.open(configFilePath);
		if (!configFileStream) 
			throw std::runtime_error("Could not open config file " + configFilePath);

	    SettingsMap items;
		std::string line;
		std::vector<std::string> tokens;
		while (std::getline(configFileStream, line)) {
			// skip empty lines and comments
			if ((line.size() <= 1) || (line[0] == '#'))
				continue;
		   
			SplitString(line, "=", tokens);

			items.insert({tokens[0], tokens[1]});
			
			tokens.clear();
		}
		configFileStream.close();

	    SettingsMap::const_iterator it;
	    SettingsMap::const_iterator end = items.end();

		// center of mass energy: REQUIRED
		// in GeV; must scale to MeV
		it = items.find("ECM");
		if (it == end)
			throw std::runtime_error("Key 'ECM' not present in configuration file.");
		ECM = std::stod(it->second);
		ECM *= 1000.0;
		S = ECM*ECM;
		LOGGER.logMessage("Using ECM=%lf", ECM);

		// name of PDF: REQUIRED
		it = items.find("PDFName");
		if (it == end)
			throw std::runtime_error("Key 'PDFName' not present in configuration file.");
		PDFName = items.at("PDFName");
		LOGGER.logMessage("Using PDF=%s", PDFName.c_str());

		// member number for PDF: NOT REQUIRED (0 default)
		it = items.find("PDFMemberNo");
		if (it == end)
			PDFMemberNo = 0;
		else
			PDFMemberNo = std::stoi(items.at("PDFMemberNo"));
		LOGGER.logMessage("\t Member ID=%d", PDFMemberNo);

		// process string: REQUIRED
		it = items.find("Process");
		if (it == end)
			throw std::runtime_error("Key 'Process' not present in configuration file.");
		Process = items.at("Process");
		LOGGER.logMessage("Process string=%s", Process.c_str());

	    it = items.find("NumXSIterations");
		if (it == end)
			throw std::runtime_error("Key 'NumXSIterations' not present in configuration file.");
		numXSIterations = std::stoi(items.at("NumXSIterations"));
		LOGGER.logMessage("Using %d iterations for cross section calculation", numXSIterations);

		
		it = items.find("AllowPhotonEmission");
		if (it == end)
			throw std::runtime_error("Key 'AllowPhotonEmission' not present in configuration file.");
		doPhotonEmission = items.at("AllowPhotonEmission").compare("yes") == 0;
		if (doPhotonEmission)
			LOGGER.logMessage("Doing photon emission\n");

		it = items.find("AllowGluonEmission");
		if (it == end)
			throw std::runtime_error("Key 'AllowGluonEmission' not present in configuration file.");
		doGluonEmission = items.at("AllowGluonEmission").compare("yes") == 0;
		if (doGluonEmission)
			LOGGER.logMessage("Doing gluon emission");

		// can't have both
		if (doPhotonEmission && doGluonEmission)
			throw std::runtime_error("Cannot do both photon and gluon emission.");

		// however at the moment we also cannot do photon emission
		// TODO: implement photon emission and remove this
		if (doPhotonEmission)
			throw std::runtime_error("Cannot do photon emission at the moment. Please do gluon emission only.");
		


	    LOGGER.logMessage("Finished loading configuration file data.");

		LHAPDF::setVerbosity(0);  // don't print anything, LHAPDF!
		pdf = LHAPDF::mkPDF(PDFName, PDFMemberNo);
	}

}; // namespace ColSim
