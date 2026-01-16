#ifndef __PARTON_SHOWER_HPP
#define __PARTON_SHOWER_HPP

#include <random>

#include "colsim/alphas.hpp"
#include "colsim/common.hpp"

namespace colsim
{
	class PartonShower
	{
	protected:
		AlphaS _alphas;

	public:
		PartonShower() : _alphas(1) {};
		~PartonShower() = default;

		struct Emission
		{
			double t, z, pT_2, m_2;
			bool generated, continue_evol;
		};


		struct Params
		{
			double Q, Qcut, r, alphas_over;
		};

	protected:
		double tgamma(double z, double alphaSOver);
		double tgamma_inv(double r, double alphaSOver);

		double zp_over(double t);
		double zm_over(double t);

		double Pqq(double z);
		double Pqq_over(double z);

		double emission_scale_func(double log_T_Q2, void* _params);
		double get_t_emission(double Q, double Qcut, double r, double tfac, double alphaSOver);

		double get_z_emission(double t, double r, double alphaSOver);
		double get_pt_2(double t, double z);
		double get_m_2(double t, double z);

		Emission generate_emission(double Q, double Qcut, double tfac, double alphaSOver);
	public:
		std::vector<Emission> evolve(double Q, double Qmin, double alphaSOver);

		AlphaS const& alphas() const { return _alphas; }
		
		


	private:
		// helper for printing the list of emissions
		// for a particle evolution
		// uses the logger's output stream
		void log_emissions_table(std::vector<Emission> const& emissions);
		
	};

	class GluonShower : public PartonShower { };
	
}; // namespace ColSim

#endif // __PARTON_SHOWER_HPP
