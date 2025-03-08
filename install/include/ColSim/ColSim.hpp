#ifndef __COLSIM_HPP
#define __COLSIM_HPP

#include "ColSim/Types.hpp"
#include "ColSim/HardProcess.hpp"
#include "ColSim/PartonShower.hpp"
#include "ColSim/Settings.hpp"

#include <memory>

namespace ColSim {

	class ColSimMain {
	private:
		std::unique_ptr<HardProcess> hardProcess;
		std::unique_ptr<PartonShower> partonShower;

		// values set after cross section calculation
		Double crossSection, crossSectionError;
		Double maxWeight;
		std::vector<Double> maxPSPoints;

		std::vector<Event> eventRecord;
		
		
	public:
		ColSimMain(const std::string& configFilePath="config.in",
				   const std::string& logFilePath="out.log");

		/** Calculates the cross section and does a few more initialization
		 *  steps in preparation for event generation.
		 */
		void start();
		
		/** Generates a single event and stores it in the event record.
		 *  Returns a const reference to this event
		 */
		const Event& generateEvent();

		/** Returns a const reference to the last event generated.
		 */
		inline const Event& getLastEvent() const {
			return eventRecord.back();
		}


		/** Finalizes event generation and consolidates the event record.
		 */
		void stop();
		
	private:		
		/** Helper function to load the corresponding hard process.
		 */
		void loadHardProcess(const std::string& processStr);

		/** Helper function to load the corresponding showering.
		 */
		void loadPartonShower(Bool doPhotonEmission, Bool doGluonEmission);
	};
	
}; // namespace ColSim;



#endif // __COLSIM_HPP
