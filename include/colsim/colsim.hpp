#ifndef __COLSIM_HPP
#define __COLSIM_HPP

#include <memory>

#include "colsim/common.hpp"
#include "colsim/hard_process.hpp"
#include "colsim/parton_shower.hpp"
#include "colsim/settings.hpp"
#include "colsim/utils.hpp"
#include "colsim/event.hpp"

namespace colsim
{
	class ColSimMain
	{
	public:
		enum InitFlag
		{
			HARD_SCATTERING,
			PARTON_SHOWERING
		};
		
	private:
		InitFlag flag;
		
		std::unique_ptr<HardProcess> _hard_process;
		std::unique_ptr<PartonShower> _parton_shower;

		// values set after cross section calculation
		double _xs, _xs_error;
		double _max_weight;
		std::vector<double> _max_ps_points;

		// the main event/emission records
		std::vector<Event> _event_record;
		std::vector<std::vector<PartonShower::Emission>> _emission_record;
		
		// list of plot points for generating plots
		 std::vector<std::vector<double>> _plot_points;
		
	public:
		ColSimMain(std::string const& log_file_path="out.log");
		~ColSimMain();

		/** Initializes either the hard scattering or parton showering process
		 *  depending on the passed flag.
		 */
	    void init(InitFlag initFlag, std::string const& config_file_path="");

		/** Calculates the cross section and does a few more initialization
		 *  steps in preparation for event generation.
		 */
		void start();
		
		/** Generates a single event and stores it in the event record.
		 */
		bool generate_event();

		/** Generates @a numEvents events and stores them in the event record.
		 */
		void generate_events(uint numEvents);

		/** Returns a const reference to the last event generated.
		 */
		inline Event& get_last_event() {
			return _event_record.back();
		}

		inline std::vector<Event> const& getEventRecord() const { return _event_record; }

		/** Takes the all of the accepted phase space points
		 *  and generates plots with them.
		 */
		void generate_plots();


		/** Finalizes event generation and consolidates the event record.
		 */
		void stop();


		/** Simple getters for the cross section and cross section errors.
		 */
		inline double cross_section() const { return _xs; }
		inline double cross_section_error() const { return _xs_error; }

		
		inline std::vector<std::vector<double>> const& plot_points() const { return _plot_points; }
		inline std::vector<std::vector<PartonShower::Emission>> const& emission_record() const { return _emission_record; }
		
	private:		
		/** Helper function to load the corresponding hard process.
		 */
		void load_hard_process(std::string const& processStr);


		// ----------------------------------------
		// Hard Process/Showering - Specific Funcs
		// ----------------------------------------
		
		/** Internal function called from @a start()
		 *  if the user passed in the @a HARD_PROCESS
		 *  initialization flag.
		 */
		void start_hard_process();

		/** Internal function called from @a start()
		 *  if the user passed in the @a PARTON_SHOWER
		 *  initialization flag.
		 */
		void start_parton_shower();

		/** Internal function called from @a generateEvent()
		 *  if the user passed in the @a HARD_PROCESS
		 *  initialization flag.
		 */
		bool generate_event_hard_process();

		/** Internal function called from @a generateEvent()
		 *  if the user passed in the @a PARTON_SHOWER
		 *  initialization flag.
		 */
		bool generate_event_parton_shower();

		/** Internal function called from @a generatePlots()
		 *  if the user passed in the @a HARD_PROCESS
		 *  initialization flag.
		 */
		void generate_plots_hard_process();

		/** Internal function called from @a generatePlots()
		 *  if the user passed in the @a PARTON_SHOWER
		 *  initialization flag.
		 */
		void generate_plots_parton_shower();

	};
	
}; // namespace ColSim;



#endif // __COLSIM_HPP
