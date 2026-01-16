#include "colsim/colsim.hpp"

#include <memory>

#include "colsim/common.hpp"
#include "colsim/utils.hpp"
#include "colsim/alphas.hpp"
#include "colsim/parton_shower.hpp"
#include "colsim/settings.hpp"
#include "colsim/phase_space.hpp"
#include "colsim/math.hpp"

namespace colsim
{
	ColSimMain::ColSimMain(std::string const& log_file_path) {
		UNUSED(log_file_path);
	    // initalize the logger with the default file path
	}

	ColSimMain::~ColSimMain() {
		log(LOG_INFO, "ColSimMain::~ColSimMain()", "Shutting down...");
	}


	void ColSimMain::init(InitFlag init_flag, std::string const& config_file_path) {
		flag = init_flag;
		
		// read in the settings from the given configuration file
		SETTINGS.load_config_file(config_file_path);

		// load the process given in the config file
		switch(init_flag) {
			case HARD_SCATTERING:
				load_hard_process(SETTINGS.process);
				break;
			case PARTON_SHOWERING:
				_parton_shower = std::unique_ptr<PartonShower>(new GluonShower());
				break;
		}
	}




	void ColSimMain::load_hard_process(std::string const& process) {
		if (process.compare("PP2Zg2ll") != 0)
			log(LOG_ERROR, "ColSimMain::load_hard_process()", "Only 'PP2Zg2ll' has been implemented.");

		_hard_process = std::unique_ptr<HardProcess>(new PP2Zg2ll());
	}


	void ColSimMain::start() {
		switch(flag) {
			case HARD_SCATTERING:
				start_hard_process();
				break;
			case PARTON_SHOWERING:
				start_parton_shower();
				break;
		}
	};
	

	bool ColSimMain::generate_event() {
		switch(flag) {
			case HARD_SCATTERING:
				return generate_event_hard_process();
			case PARTON_SHOWERING:
				return generate_event_parton_shower();
		}

		return true; // unreachable, but gcc shuts the fuck up
	}

	void ColSimMain::generate_events(uint numEvents) {
		_plot_points.clear();
		_emission_record.clear();
		while (numEvents > 0) {
			if(generate_event())
				numEvents--;
		}
	}


	void ColSimMain::generate_plots() {
		switch(flag) {
			case HARD_SCATTERING:
				generate_event_hard_process();
				break;
			case PARTON_SHOWERING:
				generate_event_parton_shower();
				break;
		}
	}

	
	void ColSimMain::stop() {
		// does nothing at the moment!
		log(LOG_INFO, "ColSimMain()::stop()", "Stopped event generation!");
	}




	bool ColSimMain::generate_event_hard_process() {
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

		// hit-or-miss to see if it actually works
		double weight_ratio = weight/_max_weight;
		double rand = rand_double();

		while (rand > weight_ratio) {
			phase_space.fill_phase_space(phase_space_points);
			res = _hard_process->dsigma(phase_space_points);
			weight = res.weight;
			for (double d : deltas)
				weight *= d;

			// ensure to calculate the new weight and another random number
			weight_ratio = weight/_max_weight;
			rand = rand_double();
		}

		// generate a list of particles
		std::vector<Particle> particles;
		_hard_process->generate_particles(particles);
		_event_record.emplace_back(Event(weight, particles));

		// plot the phase space points and the chosen additional values
		std::vector<double> plot_points(phase_space_points.begin(), phase_space_points.end());
		plot_points.append_range(res.additional_vals);
		_plot_points.emplace_back(plot_points);

		return true;
	}

	bool ColSimMain::generate_event_parton_shower() {
	    double Q0 = SETTINGS.initial_evol_e;
		double Qf = SETTINGS.evol_energy_cutoff;
		
		double scale;
		if (SETTINGS.fixed_scale)
			scale = Q0/2.0;
		else {
			scale = Qf;
		}

		AlphaS const& as = _parton_shower->alphas();
		_emission_record.emplace_back(_parton_shower->evolve(Q0, Qf, as.alphas_over(scale)));
		return true;
	}


	void ColSimMain::start_hard_process() {
		// ensure to reset ranges in the event that
		// we have changed a config file variable between calculations
		_hard_process->get_phase_space().set_ranges();
		
		log(LOG_INFO, "ColSimMain::start_hard_process()", "Calculating cross section via Monte Carlo integration...");
		HardProcessResult res = _hard_process->calculate();
		log(LOG_INFO, "ColSimMain::start_hard_process()", 
			"Result is: {:.9f} +- {:.9f} pb (picobarns)", res.result, res.error);

		_xs = res.result;
		_xs_error = res.error;
		_max_weight = res.max_weight;
		_max_ps_points = res.max_points;

		log(LOG_INFO, "ColSimMain::start_hard_process()", "Maximum weight achieved: {:.9f}", _max_weight);
	}

	void ColSimMain::start_parton_shower() {
		// nothing else to initialize
		log(LOG_INFO, "ColSimMain::start_parton_shower()", "Starting parton shower event generation.");
	}

	void ColSimMain::generate_plots_hard_process()
	{}
	void ColSimMain::generate_plots_parton_shower()
	{}

	/*
	void ColSimMain::generate_plots_hard_process() {
		log(LOG_INFO, "", "Generating plots...");

		PhaseSpace& phase_space = _hard_process->get_phase_space();

		Gnuplot plot;
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
	}
	*/

	
	/*
	void ColSimMain::generatePlots_partonShower() {
		const std::vector<std::string> plotNames{"t", "pT", "m"};

		std::vector<Double> v;
	    for (const std::vector<PartonShower::Emission>& emissions : emissionRecord) {
			for (const PartonShower::Emission& e : emissions) {
				v.push_back(std::sqrt(e.t));
				v.push_back(std::sqrt(e.pT_2));
				v.push_back(std::sqrt(e.m_2));
				plotPoints.emplace_back(v);
				v.clear();
			}
		}

		// maximum for t is directly cutoff at this value,
		// so we can fix it here
		double tMax = std::sqrt(SETTINGS.initialEvolEnergy);
		
		std::vector<Double>
			min{0.0, 0.0, 0.0},
			max{tMax, 20.0, 25.0},
			delta{tMax, 20.0, 25.0};

		
		Gnuplot plot;
		plot.setHistInfo(min, max, delta, 100);
		plot.setXLabels({"√t", "p_T", "m_{virt}"});
		plot.setYLabels({"Events", "Events", "Events"});
		plot.openDataFile("emissions.dat", plotNames);
		plot.addDataPoints(plotPoints);

		LOGGER.logMessage("Saved datafile in 'emissions.dat'");
		LOGGER.logMessage("Creating temporary Gnuplot script files for plot generation...");

		
		
		plot.plot();

		LOGGER.logMessage("Plots saved! Check the 'plots' directory to view them.");
	}
	*/
}; // namespace ColSim
