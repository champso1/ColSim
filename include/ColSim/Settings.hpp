#ifndef __SETTINGS_HPP
#define __SETTINGS_HPP

#include "ColSim/Types.hpp"

#include "LHAPDF/LHAPDF.h"

#include <unordered_map>
#include <fstream>

namespace ColSim {

	class Settings {
	public:

		// general settings
		Double ECM, S;
		std::string PDFName;
		UInt32 PDFMemberNo;
		LHAPDF::PDF* pdf;

		// hard scattering settings
		std::string Process;
		UInt32 numXSIterations;

		// parton showering settings
		Bool doPhotonEmission;
		Bool doGluonEmission;

		// values stored after the cross section
		// has been computed
		Double crossSection;
		Double crossSectionError; // standard deviation
		
	    
		Settings() : PDFName(""), Process(""), pdf(nullptr) {}
		~Settings() { delete pdf; }
		void loadConfigFile(const std::string& fileStream);

		static Settings& getInstance() {
			static Settings settings;
			return settings;
		}
	};

#define SETTINGS Settings::getInstance()
	
}; // namespace ColSim



#endif // __SETTINGS_HPP
