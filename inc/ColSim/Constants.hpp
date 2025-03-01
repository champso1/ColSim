#ifndef __CONSTANTS_HPP
#define __CONSTANTS_HPP

#include "ColSim/Types.hpp"

#include <cmath>


namespace ColSim {

#define ALPHA           7.2973525693e-3F // dimensionless
#define Z_MASS          91.1880F         // GeV/c^2
#define Z_WIDTH         2.4414F          // GeV
#define MAGIC_FACTOR    3.893793721e8F   // GeV^2 pb
#define FERMI_CONSTANT  1.1663788e-5F    // GeV^-2
#define WEINBERG_ANGLE  0.222246F        // dimensionless 
#define PI              3.141592653F     // retroencabulators / football field


	constexpr const Double Z_MASS_2 = Z_MASS*Z_MASS;
	constexpr const Double Z_WIDTH_2 = Z_WIDTH*Z_WIDTH;


	// some more process-specific constants
#define Q_MIN           60.0F            // minimum energy cutoff
#define MASS_TR         60.0F            // transformed mass
#define WIDTH_TR        60.0F            // transformed width



	constexpr const double Q_MIN_2 = Q_MIN * Q_MIN;
	constexpr const double MASS_TR_2 = MASS_TR * MASS_TR;
	constexpr const double WIDTH_TR_2 = WIDTH_TR * WIDTH_TR;
	
};





#endif // __CONSTANTS_HPP
