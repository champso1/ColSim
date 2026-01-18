#ifndef __SETTINGS_HPP
#define __SETTINGS_HPP

#include <unordered_map>
#include <memory>

#include "LHAPDF/LHAPDF.h"

namespace colsim
{
	inline auto lhapdf_pdf_deleter = [](LHAPDF::PDF* pdf) { delete pdf; };
	using lhapdf_pdf_deleter_type = decltype(lhapdf_pdf_deleter);
	using lhapdf_pdf_type = std::unique_ptr<LHAPDF::PDF, lhapdf_pdf_deleter_type>;

	struct Settings final
	{
		using value_type = std::unordered_map<std::string, std::string>;

		// general settings
		double ecm, s;
		std::string pdf_name{};
		int pdf_mem;
	    lhapdf_pdf_type pdf{nullptr, lhapdf_pdf_deleter};

		// hard scattering settings
		std::string process{};
		int num_iterations;
		double min_cutoff_energy, min_cutoff_energy_2;
		double trans_energy, trans_energy_2;

		// parton showering settings
		double initial_evol_e, initial_evol_e_2;
		bool fixed_scale;
		double evol_energy_cutoff;
		
		// values stored after the cross section
		// has been computed
		double xs;
		double xs_error;

		value_type settings{};
		
		Settings() = default;
		~Settings() = default;

		void load_config_file(std::string const& config_file_path);
		void read_string(std::string const& str);

		static Settings& instance()
		{
			static Settings settings;
			return settings;
		}

	private:
		bool does_key_exist(value_type::const_iterator& it, std::string const& key);
	};

#define SETTINGS Settings::instance()
	
}; // namespace ColSim



#endif // __SETTINGS_HPP
