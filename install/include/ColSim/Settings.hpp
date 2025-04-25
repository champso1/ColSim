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
		Double minCutoffEnergy, minCutoffEnergy_2;
		Double transEnergy, transEnergy_2;

		// parton showering settings
		Double initialEvolEnergy, initialEvolEnergy_2;
		Bool fixedScale;
		Double evolEnergyCutoff;
		

		// values stored after the cross section
		// has been computed
		Double crossSection;
		Double crossSectionError; // standard deviation

		
		typedef std::unordered_map<std::string, std::string> SettingsMap;
		SettingsMap settings;
		
	    
		Settings() : PDFName(""), pdf(nullptr), Process("") {}
		~Settings() { delete pdf; }
		void loadConfigFile(const std::string& fileStream);

		/** Reads @a str as if it were a line in the configuration file.
		 *  Can be used to set values after the configuration file,
		 *  like if multiple runs are to be done during one call of the executable.
		 */
		void readString(const std::string& str);

		static Settings& getInstance() {
			static Settings settings;
			return settings;
		}

	private:
		/** Helper for determining whether a key is inside the configuration file.
		 *  If found, returns true and sets @a it to point to the item.
		 */
		Bool doesKeyExist(SettingsMap::const_iterator& it, const std::string& key);
	};

#define SETTINGS Settings::getInstance()
	
}; // namespace ColSim



#endif // __SETTINGS_HPP
