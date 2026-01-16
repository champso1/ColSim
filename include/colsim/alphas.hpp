#ifndef __CONSTANTS_HPP
#define __CONSTANTS_HPP

#include "colsim/common.hpp"
#include "colsim/settings.hpp"

namespace colsim
{
	/** Class for determining values of the QCD running coupling
	 */
	class AlphaS
	{
		uint _order; // perturbative order
		double _alphas_b, _alphas_c; // calced values of alphaS at B and C masses
		double _fixed_scale_val;

	public:
		static constexpr double ALPHAS_Z = 0.118;

		AlphaS(uint _order)
			: _order{_order},
			  _alphas_b{calc_alphas(B_MASS_2)}, _alphas_c{calc_alphas(C_MASS_2)},
			  _fixed_scale_val{SETTINGS.initial_evol_e/2.0}
		{}

	    double beta0(uint nf) const;
		double beta1(uint nf) const;
		double calc_alphas(double t) const;
		double alphas_scale(double t, double z) const;
		double alphas_actual(double t, double z) const;
		double alphas_over(double scale) const;

	private:
		double calc_alphas_LO(double t) const;
		double calc_alphas_NLO(double t) const;
	};
};





#endif // __CONSTANTS_HPP
