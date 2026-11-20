#ifndef __COLSIM_HPP
#define __COLSIM_HPP

#include <memory>

#include "colsim/common.hpp"
#include "colsim/hard_process.hpp"
#include "colsim/parton_shower.hpp"
#include "colsim/event.hpp"

namespace colsim
{
	class ColSimMain final
	{
	public:
		enum class InitFlag
		{
			HARD_SCATTERING,
			PARTON_SHOWERING
		};

		struct CrossSectionInfo final
		{
			double xs, xs_error;
			double max_weight;
			std::vector<double> _max_ps_points;
		};
		
	private:
		InitFlag _init_flag;
		
		std::unique_ptr<HardProcess> _hard_process{};

		// values set after cross section calculation
	    CrossSectionInfo _xs_info{};

		// the main event/emission records
		std::vector<Event> _event_record{};

	public:
		ColSimMain(InitFlag initFlag);

		/** Calculates the cross section and does a few more initialization
		 *  steps in preparation for event generation.
		 */
		void start();
		
		/** Generates a single event and stores it in the event record.
		 *  @returns reference to generated event
		 */
		Event const& generate_event();

		/** Generates @a num_events events and stores them in the event record.
		 *  @returns reference to stored event record
		 */
		std::vector<Event> const& generate_events(uint num_events);

		/** Simple getters for the cross section and cross section errors.
		 */
		inline CrossSectionInfo const& cross_section_info() const { return _xs_info; }

		/**
		 *  @brief spits out 0 plots right now lolol
		 */
		void generate_plots();
		
	private:		
		// ----------------------------------------
		// Hard Process/Showering - Specific Funcs
		// ----------------------------------------
		
		/** Internal function called from @a start()
		 *  if the user passed in the @a HARD_PROCESS
		 *  initialization flag.
		 */
		void start_hard_process();

		/** Internal function called from @a generateEvent()
		 *  if the user passed in the @a HARD_PROCESS
		 *  initialization flag.
		 */
		bool generate_event_hard_process();


		/** Internal function called from @a generatePlots()
		 *  if the user passed in the @a HARD_PROCESS
		 *  initialization flag.
		 */
		void generate_plots_hard_process();

	};
	
}; // namespace ColSim;



#endif // __COLSIM_HPP
