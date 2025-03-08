#include "ColSim/Constants.hpp"
#include "ColSim/Settings.hpp"


namespace ColSim {

	Double AlphaS::Beta0(UInt8 nf) {
		Double Nf = static_cast<Double>(nf);
		return (11.0/6.0)*CA - (2.0/3.0)*TR*Nf;
	}

	Double AlphaS::Beta1(UInt8 nf) {
		Double Nf = static_cast<Double>(nf);
		return (17.0/6.0)*std::pow(CA,2) - ((5.0/3.0)*CA + CF)*TR*Nf;
	}

	Double AlphaS::calcAlphaS_LO(Double t) {
		// determine reference scale given t
		Double tref, asref;
		Double beta0;

		if (t >= B_MASS_2) {
			// if above B mass, use Z mass
			tref = Z_MASS_2;
			asref = ALPHAS_Z;
			// 5 flavors
			beta0 = Beta0(5) / (2.0*PI);
		} else if (t >= C_MASS_2) {
			// use B
			tref = B_MASS_2;
			asref = alphaS_B;
			// 4 flavors now
			beta0 = Beta0(4) / (2.0*PI);
		} else {
			// otherwise just use C
			tref = C_MASS_2;
			asref = alphaS_C;
			// 3 flavors
			beta0 = Beta0(3) / (2.0*PI);
		}

		return 1.0/((1.0/asref) + beta0*std::log(t/tref));
	}

	Double AlphaS::calcAlphaS_NLO(Double t) {
		// determine reference scale given t
		Double tref, asref;
		Double beta0, beta1;

		if (t >= B_MASS_2) {
			// if above B mass, use Z mass
			tref = Z_MASS_2;
			asref = ALPHAS_Z;
			// 5 flavors
			beta0 = Beta0(5) / (2.0*PI);
			beta1 = Beta1(5) / std::pow(2.0*PI,2);
		} else if (t >= C_MASS_2) {
			// use B
			tref = B_MASS_2;
			asref = alphaS_B;
			// 4 flavors now
			beta0 = Beta0(4) / (2.0*PI);
			beta1 = Beta1(4) / std::pow(2.0*PI,2);
		} else {
			// otherwise just use C
			tref = C_MASS_2;
			asref = alphaS_C;
			// 3 flavors
			beta0 = Beta0(3) / (2.0*PI);
			beta1 = Beta1(3) / std::pow(2.0*PI,2);
		}

		Double term = 1.0 + beta0*asref*std::log(t/tref);
		return (asref/term) * (1.0 - (beta1/beta0)*asref*std::log(term)/term);
	}


	Double AlphaS::calcAlphaS(Double t) {
		switch (order) {
			case 0: return calcAlphaS_LO(t);
			case 1: return calcAlphaS_NLO(t);
		}
		return 0.0; // just in case
	}


	Double AlphaS::alphaSActual() {
	    const LHAPDF::PDF* pdf = SETTINGS.pdf;
		return pdf->alphasQ(Q_SCALE_FIXED) / (2.0*PI);
	}

	Double AlphaS::alphaSOver() {
		return alphaSActual();
	}

	
}; // namespace ColSim
