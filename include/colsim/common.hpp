#ifndef __COMMON_HPP
#define __COMMON_HPP

#include <filesystem>

#define COLSIM_JOIN0(a, b) a ## b
#define COLSIM_JOIN(a, b) COLSIM_JOIN0(a, b)

#ifndef UNUSED
#define UNUSED(x) ((void)x)
#endif

namespace colsim
{

    namespace fs
    {
        using namespace std::filesystem;
    }

    using uint = unsigned;

    // group theoretical constants
    constexpr double NC = 3.0;
    constexpr double TR = (1.0/2.0);
    constexpr double CA = NC;
	constexpr double CF = ((NC*NC) - 1.0) / (2.0*NC);

// other actual constant values (non-configurable)
    constexpr double ALPHA           = 7.2973525693e-3; // QED/EW coupling constant
    constexpr double Z_MASS          = 91.1880;         // mass of Z boson
    constexpr double Z_WIDTH         = 2.4414;          // decay width of Z boson
    constexpr double C_MASS          = 1.273;           // c quark mass
    constexpr double B_MASS          = 4.183;           // b quark mass
    constexpr double MAGIC_FACTOR    = 3.893793721e8;   // conversion factor
    constexpr double FERMI_CONSTANT  = 1.1663788e-5;
    constexpr double WEINBERG_ANGLE  = 0.222246;
    constexpr double PI              = 3.141592653;
	
	// squares of some of the above values
	constexpr double Z_MASS_2 = Z_MASS*Z_MASS;
    constexpr double Z_WIDTH_2 = Z_WIDTH * Z_WIDTH;
	constexpr double C_MASS_2 = C_MASS*C_MASS;
	constexpr double B_MASS_2 = B_MASS*B_MASS;
}

#endif // __COMMON_HPP