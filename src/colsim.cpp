#include "colsim/colsim.hpp"

#include <memory>

#include "colsim/common.hpp"
#include "colsim/alphas.hpp"
#include "colsim/phase_space.hpp"
#include "colsim/math.hpp"
#include "colsim/logger.hpp"
#include "colsim/event.hpp"

namespace colsim
{
	ColSimMain::ColSimMain(InitFlag init_flag)
		: _init_flag{init_flag}
	{
		_hard_process = std::make_unique<PP2Zg2ll>();
	}

	void ColSimMain::start()
	{
		start_hard_process();
	};
	

	Event const& ColSimMain::generate_event()
	{
		generate_event_hard_process();
		return _event_record.back();
	}

	std::vector<Event> const& ColSimMain::generate_events(uint num_events)
	{
		while (num_events > 0) {
			if(generate_event().valid())
				num_events--;
		}
		return _event_record;
	}
	

	void ColSimMain::start_hard_process()
	{
		// ensure to reset ranges in the event that
		// we have changed a config file variable between calculations
		_hard_process->get_phase_space().set_ranges();
		
		log(LOG_INFO, "ColSimMain::start_hard_process()", "Calculating cross section via Monte Carlo integration...");
		HardProcessResult res = _hard_process->calculate();
		log(LOG_INFO, "ColSimMain::start_hard_process()", 
			"Result is: {:.9f} +- {:.9f} pb (picobarns)", res.result, res.error);

		_xs_info = CrossSectionInfo{
			.xs = res.result,
			.xs_error = res.error,
			.max_weight = res.max_weight,
			.max_ps_points = res.max_points
		};

		log(LOG_INFO, "ColSimMain::start_hard_process()", "Maximum weight achieved: {:.9f}", _xs_info.max_weight);
	}

	
	void ColSimMain::generate_event_hard_process()
	{
		PhaseSpace phase_space = _hard_process->get_phase_space();

		std::vector<double> phase_space_points;
		phase_space.fill_phase_space(phase_space_points);
		std::vector<double> const& deltas = phase_space.deltas();
		
		HardProcess::Result res = _hard_process->dsigma(phase_space_points);
		if (res == HardProcess::Result::invalid_result())
			return false;
		
		double weight = res.weight;
		// scale the weight
		for (double d : deltas)
			weight *= d;

		double max_weight = _xs_info.max_weight;

		// hit-or-miss to see if it actually works
		double weight_ratio = weight/max_weight;
		double rand = rand_double();

		while (rand > weight_ratio) {
			phase_space.fill_phase_space(phase_space_points);
			res = _hard_process->dsigma(phase_space_points);
			weight = res.weight;
			for (double d : deltas)
				weight *= d;

			// ensure to calculate the new weight and another random number
			weight_ratio = weight/max_weight;
			rand = rand_double();
		}

		// generate a list of particles
		std::vector<Particle> particles;
		_hard_process->generate_particles(particles);
		_event_record.emplace_back(Event(weight, particles));
	}

	void ColSimMain::generate_plots()
	{
	    generate_plots_hard_process();
	}
	
	void ColSimMain::generate_plots_hard_process()
	{
		log(LOG_WARNING, "ColSimMain::generate_plots_hard_process()", "Not implemented this sheisse yet.");
		
		/*
		log(LOG_INFO, "", "Generating plots...");

		PhaseSpace& phase_space = _hard_process->get_phase_space();

		Plotter<double> plot;
		
		Plotter plot;
		plot.setHistInfo(phase_space.getMins(),
						 phase_space.getMaxes(),
						 phase_space.getDeltas(), 100);
		std::vector<std::string> plotColNames(phase_space.getNames());
		plot.openDataFile("events.dat", plotColNames);
	    plot.addDataPoints(plotPoints);

		LOGGER.logMessage("Saved datafile in 'events.dat'");
		LOGGER.logMessage("Creating temporary Gnuplot script files for plot generation...");

		// plot.setTitles(phase_space.getTitles());
		plot.setXLabels(phase_space.getXLabels());
		plot.setYLabels(phase_space.getYLabels());
		
		plot.plot();

		LOGGER.logMessage("Plots saved! Check the 'plots' directory to view them.");
		*/
	}
	
}; // namespace ColSim
