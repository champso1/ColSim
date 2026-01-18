#include "colsim/settings.hpp"
#include "colsim/utils.hpp"

#include <cmath>
#include <ranges>
#include <string>
#include <vector>
#include <fstream>
#include <unordered_map>

#include "LHAPDF/LHAPDF.h"

namespace colsim
{

	bool Settings::does_key_exist(Settings::value_type::const_iterator& it, std::string const& key) {
		bool found;
		it = settings.find(key);
		found = (it != settings.end());

		return found;
	}
	

	void Settings::load_config_file(std::string const& config_file_path) {
		if (config_file_path.length() > 0) {
			std::ifstream config_file_stream;
			config_file_stream.open(config_file_path);
			if (!config_file_stream) 
				log(LOG_ERROR, "Settings::load_config_file", "Could not open config file '{}'", config_file_path);
				
			std::string line;
			std::vector<std::string> tokens;
			while (std::getline(config_file_stream, line)) {
				// skip empty lines and comments
				if ((line.size() <= 1) || (line[0] == '#'))
					continue;
				
				tokens = 
					line
					| std::ranges::views::split('=')
					| std::ranges::to<std::vector<std::string>>();

				settings.insert({tokens[0], tokens[1]});
			
				tokens.clear();
			}
			config_file_stream.close();
		}

	    value_type::const_iterator it;

		// center of mass energy
		// in TeV; must scale to GeV
		if (does_key_exist(it, "ECM")) {
			ecm = std::stod(it->second);
			if (ecm < 5.0 || ecm > 20.0)
				log(LOG_WARNING, "Settings::load_config_file()", "ECM value of {} is abnormal.", ecm);
			ecm *= 1000.0;
			s = ecm*ecm;
		} else {
			ecm = 14000.0;
			s = ecm*ecm;
		}
		log(LOG_INFO, "Settings::load_config_file()", "Using ECM={}", ecm);

		// name of PDF
		if(does_key_exist(it, "PDFName"))
			pdf_name = it->second;
		else
			pdf_name = "CT18NNLO";
		LHAPDF::setVerbosity(0);
		pdf = std::unique_ptr<LHAPDF::PDF, lhapdf_pdf_deleter_type>(LHAPDF::mkPDF(pdf_name, 0), lhapdf_pdf_deleter);

		// process string
		process = "PP2Zg2ll";
		log(LOG_INFO, "Settings::load_config_file()", "Process string={}", process);

		// number if iterations for cross section calculation: REQUIRED
		if(does_key_exist(it, "NumXSIterations"))
			num_iterations = std::stoi(it->second);
		else
			num_iterations = 1000000;
			
		log(LOG_INFO, "Settings::load_config_file()", "Using {} iterations for cross section calculation", num_iterations);

		// cutoff energy for phase space generation
		if(does_key_exist(it, "MinCutoffEnergy")) {
			min_cutoff_energy = std::stod(it->second);
			if (min_cutoff_energy < 1.0)
				log(LOG_ERROR, "Settings::load_config_file()", "The specified cutoff energy {:.2f} is less than the absolute minimum of 1.0.", min_cutoff_energy);
			else if (min_cutoff_energy > 500.0)
				log(LOG_WARNING, "Settings::load_config_file()", "The specified cutoff energy {:.2f} is higher than ordinary.", min_cutoff_energy);
			min_cutoff_energy_2 = min_cutoff_energy*min_cutoff_energy;
		} else {
			min_cutoff_energy = 60.0;
			min_cutoff_energy_2 = 60.0*60.0;
		}
			log(LOG_INFO, "Settings::load_config_file()", "Setting minimum cutoff energy for cross section calculation to {}", min_cutoff_energy);

		// transformation mass/width (also for phase space points): REQUIRED
		if(does_key_exist(it, "TransformationEnergy")) {
			trans_energy = std::stod(it->second);
			if (trans_energy < min_cutoff_energy) {
				log(LOG_ERROR, "Settings::load_config_file()", "The specified transformation energy {} cannot be less than Q_min.", trans_energy);
			} else if (min_cutoff_energy > std::sqrt(ecm))
				log(LOG_ERROR, "Settings::load_config_file()", "The specified transformation energy {} cannot be larger than sqrt(ECM).", trans_energy);
			trans_energy_2 = trans_energy*trans_energy;
		} else {
			trans_energy = min_cutoff_energy;
			trans_energy_2 = trans_energy*trans_energy;
		}
		log(LOG_INFO, "Settings::load_config_file()", "Setting transformation mass/energy to {}", trans_energy);

		// initial evolution scale for parton showering: REQUIRED
		if(does_key_exist(it, "InitialEvolEnergy")) {
			initial_evol_e = std::stod(it->second);
			if (initial_evol_e < 100.0 || initial_evol_e > 5000.0)
				log(LOG_WARNING, "Settings::load_config_file()", "Initial evolution energy of {} is out of ordinary range of [100.0,5000.0] GeV.", initial_evol_e);
			initial_evol_e_2 = initial_evol_e*initial_evol_e;
		} else {
			initial_evol_e = 1000.0;
			initial_evol_e_2 = 1000.0*1000.0;
		}
		log(LOG_INFO, "Settings::load_config_file()", "Will start parton evolution at {}", initial_evol_e);


		// used fixed renormalization scale
		if(does_key_exist(it, "FixedScale")) {
			fixed_scale = settings.at("FixedScale").compare("Yes") == 0;
		} else {
			fixed_scale = true;
		}

		if (fixed_scale) {
			log(LOG_INFO, "Settings::load_config_file()", "Using fixed scale (mass of Z boson) for parton evolution.");
			log(LOG_WARNING, "Settings::load_config_file()", "A fixed scale misses out on some higher order effects.");
		}
		else {
			log(LOG_INFO, "Settings::load_config_file()", "Using variable scale for parton evolution.");
		}

		if(does_key_exist(it, "EvolutionEnergyCutoff")) {
			evol_energy_cutoff = std::stod(it->second);
			if (evol_energy_cutoff < 1.0) 
				log(LOG_ERROR, "Settings::load_config_file()", "Minimum cutoff energy for parton evolution MUST be greater than 1.0.");
			if (evol_energy_cutoff >= 100.0)
				log(LOG_WARNING, "Settings::load_config_file()", "Specified cutoff energy for parton evolution is abnormally high.");
		} else {
			evol_energy_cutoff = 1.0;
		}
		log(LOG_INFO, "Settings::load_config_file()", "Will cut off parton evolution at %.3lf", evol_energy_cutoff);


		log(LOG_INFO, "Settings::load_config_file()", "-----------------------------------------");
	    log(LOG_INFO, "Settings::load_config_file()", "Finished loading configuration file data.");
		log(LOG_INFO, "Settings::load_config_file()", "-----------------------------------------");
	}

}; // namespace ColSim
