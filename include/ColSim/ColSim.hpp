#ifndef __COLSIM_HPP
#define __COLSIM_HPP

#include "ColSim/Types.hpp"
#include "ColSim/HardProcess.hpp"
#include "ColSim/PartonShower.hpp"
#include "ColSim/Settings.hpp"
#include "ColSim/Logger.hpp"

#include <memory>

namespace ColSim {

	class ColSimMain {
	public:
		enum InitFlag {
			HARD_SCATTERING,
			PARTON_SHOWERING
		};
		
	private:
		InitFlag flag;
		
		std::unique_ptr<HardProcess> hardProcess;
		std::unique_ptr<PartonShower> partonShower;

		// values set after cross section calculation
		Double crossSection, crossSectionError;
		Double maxWeight;
		std::vector<Double> maxPSPoints;

		// the main event/emission records
		std::vector<Event> eventRecord;
		std::vector<std::vector<PartonShower::Emission>> emissionRecord;
		
		// list of plot points for generating plots
		 std::vector<std::vector<Double>> plotPoints;
		
	public:
		ColSimMain(const std::string& logFilePath="out.log");
		~ColSimMain();

		

		/** Initializes either the hard scattering or parton showering process
		 *  depending on the passed flag.
		 */
	    void init(InitFlag initFlag, const std::string& configFilePath="");

		/** Calculates the cross section and does a few more initialization
		 *  steps in preparation for event generation.
		 */
		void start();
		
		/** Generates a single event and stores it in the event record.
		 */
		Bool generateEvent();

		/** Generates @a numEvents events and stores them in the event record.
		 */
		void generateEvents(UInt32 numEvents);

		/** Returns a const reference to the last event generated.
		 */
		inline Event& getLastEvent() {
			return eventRecord.back();
		}

		inline const std::vector<Event>& getEventRecord() const { return eventRecord; }

		/** Takes the all of the accepted phase space points
		 *  and generates plots with them.
		 */
		void generatePlots();


		/** Finalizes event generation and consolidates the event record.
		 */
		void stop();


		/** Simple getters for the cross section and cross section errors.
		 */
		inline Double getCrossSection() const { return crossSection; }
		inline Double getCrossSectionError() const { return crossSectionError; }

		
		inline const
		std::vector<std::vector<Double>>&
		getPlotPoints() const {
			return plotPoints;
		}

		inline const
		std::vector<std::vector<PartonShower::Emission>>&
		getEmissionRecord() const {
			return emissionRecord;
		}
		
	private:		
		/** Helper function to load the corresponding hard process.
		 */
		void loadHardProcess(const std::string& processStr);


		// ----------------------------------------
		// Hard Process/Showering - Specific Funcs
		// ----------------------------------------
		
		/** Internal function called from @a start()
		 *  if the user passed in the @a HARD_PROCESS
		 *  initialization flag.
		 */
		void start_hardProcess();

		/** Internal function called from @a start()
		 *  if the user passed in the @a PARTON_SHOWER
		 *  initialization flag.
		 */
		void start_partonShower();

		/** Internal function called from @a generateEvent()
		 *  if the user passed in the @a HARD_PROCESS
		 *  initialization flag.
		 */
		Bool generateEvent_hardProcess();

		/** Internal function called from @a generateEvent()
		 *  if the user passed in the @a PARTON_SHOWER
		 *  initialization flag.
		 */
		Bool generateEvent_partonShower();

		/** Internal function called from @a generatePlots()
		 *  if the user passed in the @a HARD_PROCESS
		 *  initialization flag.
		 */
		void generatePlots_hardProcess();

		/** Internal function called from @a generatePlots()
		 *  if the user passed in the @a PARTON_SHOWER
		 *  initialization flag.
		 */
		void generatePlots_partonShower();

	};
	
}; // namespace ColSim;



#endif // __COLSIM_HPP
