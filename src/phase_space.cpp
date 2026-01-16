#include "colsim/phase_space.hpp"

#include <random>

#include "colsim/settings.hpp"
#include "colsim/math.hpp"

namespace colsim
{
	void PhaseSpace::fill_delta()
	{
		_delta.clear();
		for (uint i=0; i<_min.size(); i++)
			_delta.push_back(_max.at(i) - _min.at(i));
	}

	void PhaseSpace::fill_phase_space(std::vector<double>& vec)
	{
		vec.clear();
	    for (uint i=0; i<_num_dims; i++) {
			double val = rand_double()*_delta[i] + _min[i];
			vec.push_back(val);
		}
	}

	
	
	
	PhaseSpace_TauYCosth::PhaseSpace_TauYCosth()
		: PhaseSpace(
			3,
			{"cos(theta)", "rho", "y", "Q", "x1"},
			{"cos(θ)", "ρ", "rapidity", "Q", "x_1"},
			{"cos(θ)", "ρ", "y", "Q [GeV]", "x_1"},
			{"Events", "Events", "Events", "Events", "Events"})
	{
		
	    set_ranges();
	}

	void PhaseSpace_TauYCosth::set_ranges()
	{
		const double S = SETTINGS.s;
		const double E_TR_2 = SETTINGS.trans_energy_2;
		const double Q_MIN_2 = SETTINGS.min_cutoff_energy_2;

		_rho_min = std::atan(Q_MIN_2-E_TR_2) / E_TR_2;
		_rho_max = std::atan((S-E_TR_2) / (E_TR_2));
		_drho = _rho_max - _rho_min;

		_min.clear();
		_max.clear();
		_min.push_back(-1.0); _max.push_back(1.0);      // cosTheta
		_min.push_back(_rho_min); _max.push_back(_rho_max); // rho
		_min.push_back(0.0); _max.push_back(1.0);       // y

		// Q and x1, which aren't themselves
		// independent phase space variables
		// but are variables of interest later on
		_min.push_back(60.0); _max.push_back(500.0);  // Q
		_min.push_back(0.0); _max.push_back(1.0);    // x

		fill_delta();
	}


	PhaseSpace_EtEta::PhaseSpace_EtEta()
		: PhaseSpace(2, {"E_t", "eta"})
	{
		set_ranges();
	}


	void PhaseSpace_EtEta::set_ranges()
	{
		_min.push_back(0.1); _max.push_back(500.0); // E_t
		_min.push_back(-5.0); _max.push_back(5.0);    // eta

		fill_delta();
	}
}; // namespace ColSim
