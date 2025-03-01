#ifndef __COLSIM_HPP
#define __COLSIM_HPP

#include "ColSim/Types.hpp"
#include "ColSim/HardProcess.hpp"

#include <memory>

namespace ColSim {

	class ColSimMain {
	private:
		std::unique_ptr<HardProcess> hardProcess;
	public:
		ColSimMain(const std::string& configFilePath="./res/config.in");

		void CalcCrossSection(Double& result, Double& error);
		void GenerateEvent();
		
	private:
		// parse the hard process string given in the configuration file
		void loadHardProcess(const std::string& processStr);
	};
	
}; // namespace ColSim;



#endif // __COLSIM_HPP
