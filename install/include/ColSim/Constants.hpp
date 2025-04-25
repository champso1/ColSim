#ifndef __CONSTANTS_HPP
#define __CONSTANTS_HPP

#include "ColSim/Types.hpp"
#include "ColSim/Settings.hpp"

#include <cmath>


namespace ColSim {

// group theoretical constants
#define NC 3.0
#define TR (1.0/2.0)
#define CA NC
	constexpr const Double CF = (((NC*NC) - 1.0) / (2.0*NC));

// other actual constant values (non-configurable)
#define ALPHA           7.2973525693e-3F // QED/EW coupling constant
#define Z_MASS          91.1880F         // mass of Z boson
#define Z_WIDTH         2.4414F          // decay width of Z boson
#define C_MASS          1.273F           // c quark mass
#define B_MASS          4.183F           // b quark mass
#define MAGIC_FACTOR    3.893793721e8F   // conversion factor
#define FERMI_CONSTANT  1.1663788e-5F
#define WEINBERG_ANGLE  0.222246F
#define PI              3.141592653F
	
	// squares of some of the above values
	constexpr const Double Z_MASS_2 = Z_MASS*Z_MASS;
    constexpr const Double Z_WIDTH_2 = Z_WIDTH * Z_WIDTH;
	constexpr const Double C_MASS_2 = C_MASS*C_MASS;
	constexpr const Double B_MASS_2 = B_MASS*B_MASS;

	
#define ALPHAS_Z        0.118F // known value of alpha_s at z mass (non-configurable)

	// hard scattering
	// #define Q_MIN           60.0F            // minimum energy cutoff
	// #define MASS_TR         60.0F            // transformed mass
	// #define WIDTH_TR        60.0F            // transformed width
	
	// constexpr const Double Q_MIN_2 = Q_MIN * Q_MIN;
	// constexpr const Double MASS_TR_2 = MASS_TR * MASS_TR;
	// constexpr const Double WIDTH_TR_2 = WIDTH_TR * WIDTH_TR;

	
	// showering
	// #define Q_CUTOFF        1.0F             // minimum/final energy
	// #define Q_HARD          1000.0F          // max/initial energy
	// #define N_EVOLUTIONS    1000
	
	// constexpr const Double Q_SCALE_FIXED = Q_HARD/2.0;


	/** Class for determining values of the QCD running coupling
	 */
	class AlphaS {
	private:
		UInt8 order; // perturbative order
		Double alphaS_B, alphaS_C; // calced values of alphaS at B and C masses
		Double fixedScaleVal;

	public:
		AlphaS(UInt8 _order)
			: order(_order),
			  alphaS_B(calcAlphaS(B_MASS_2)), alphaS_C(calcAlphaS(C_MASS_2)),
			  fixedScaleVal(SETTINGS.initialEvolEnergy/2.0)
		{}

		/** Determines the first beta coefficient (i.e. at LO)
		 *  given the number of current massive/active flavors
		 */
	    Double Beta0(UInt8 nf) const;
		
		/** Determines the second beta coefficient (i.e. at NLO)
		 *  given the number of current massive/active flavors
		 */
		Double Beta1(UInt8 nf) const;

		/** Calculates the value of the QCD running coupling
		 *  at the given energy scale t either at LO or NLO
		 *  given by the initialization parameter @a order
		 */
		Double calcAlphaS(Double t) const;


		/** Determines the current scale of alphaS
		 */
		Double alphaSScale(Double t, Double z) const;

		/** Returns the actual value of the QCD running coupling
		 *  from the current PDF set
		 */
		Double getAlphaSActual(Double t, Double z) const;

		/** Returns an overestimate for the value of the QCD running coupling
		 *  to be used in simplifying and making the parton showering
		 *  algorithm more efficient.
		 */
		Double getAlphaSOver(Double scale) const;



	private:
		/** Internal function that calculates the QCD running coupling
		 *  at leading order
		 */
		Double calcAlphaS_LO(Double t) const;

		/** Internal function that calculates the QCD running coupling
		 *  at next-to leading order
		 */
		Double calcAlphaS_NLO(Double t) const;
	};


	
};





#endif // __CONSTANTS_HPP
