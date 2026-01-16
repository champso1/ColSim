#include "colsim/alphas.hpp"

namespace colsim
{
    
    double AlphaS::beta0(uint nf) const {
		
        double Nf = static_cast<
        double>(nf);
		return (11.0/6.0)*CA - (2.0/3.0)*TR*Nf;
	}

	
    double AlphaS::beta1(uint nf) const {
		
        double Nf = static_cast<
        double>(nf);
		return (17.0/6.0)*std::pow(CA,2) - ((5.0/3.0)*CA + CF)*TR*Nf;
	}

	
    double AlphaS::calc_alphas_LO(double t) const {
		// determine reference scale given t
		
        double tref, asref;
        double b0;

		if (t >= B_MASS_2) {
			// if above B mass, use Z mass
			tref = Z_MASS_2;
			asref = ALPHAS_Z;
			// 5 flavors
			b0 = beta0(5) / (2.0*PI);
		} else if (t >= C_MASS_2) {
			// use B
			tref = B_MASS_2;
			asref = _alphas_b;
			// 4 flavors now
			b0 = beta0(4) / (2.0*PI);
		} else {
			// otherwise just use C
			tref = C_MASS_2;
			asref = _alphas_c;
			// 3 flavors
			b0 = beta0(3) / (2.0*PI);
		}

		return 1.0/((1.0/asref) + b0*std::log(t/tref));
	}

	
    double AlphaS::calc_alphas_NLO( double t) const {
		// determine reference scale given t
		
        double tref, asref;
        double b0, b1;

		if (t >= B_MASS_2) {
			// if above B mass, use Z mass
			tref = Z_MASS_2;
			asref = ALPHAS_Z;
			// 5 flavors
			b0 = beta0(5) / (2.0*PI);
			b1 = beta1(5) / std::pow(2.0*PI,2);
		} else if (t >= C_MASS_2) {
			// use B
			tref = B_MASS_2;
			asref = _alphas_b;
			// 4 flavors now
			b0 = beta0(4) / (2.0*PI);
			b1 = beta1(4) / std::pow(2.0*PI,2);
		} else {
			// otherwise just use C
			tref = C_MASS_2;
			asref = _alphas_c;
			// 3 flavors
			b0 = beta0(3) / (2.0*PI);
			b1 = beta1(3) / std::pow(2.0*PI,2);
		}

		
        double term = 1.0 + b0*asref*std::log(t/tref);
		return (asref/term) * (1.0 - (b1/b0)*asref*std::log(term)/term);
	}


	double AlphaS::alphas_scale(double t, double z) const {
		if (SETTINGS.fixed_scale)
			return _fixed_scale_val;
		else
			return z*(1.0-z)*sqrt(t);

		
	}


	double AlphaS::calc_alphas(double t) const {
		switch (_order) {
			case 0: return calc_alphas_LO(t);
			case 1: return calc_alphas_NLO(t);
		}
		return 0.0; // just avoid error
	}


	double AlphaS::alphas_actual(double t, double z) const {
		double scale = alphas_scale(t, z);
		if (scale < SETTINGS.evol_energy_cutoff)
			scale = SETTINGS.evol_energy_cutoff;
		return calc_alphas(t)/(2.0*M_PI);
		
	}

	double AlphaS::alphas_over(double t) const {
		double scale;
		if (SETTINGS.fixed_scale)
			scale = _fixed_scale_val;
		else
			scale = t;
		return calc_alphas(scale);
	}
}