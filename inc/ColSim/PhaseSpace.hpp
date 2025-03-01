#ifndef __PHASE_SPACE_HPP
#define __PHASE_SPACE_HPP

#include "ColSim/Types.hpp"
#include "ColSim/Constants.hpp"

#include <array>
#include <random>

namespace ColSim {

	class PhaseSpace {
	protected:
		Double * min;
		Double * max;
		Double * delta;
		UInt32 numDims;
		
		std::default_random_engine gen;
		std::uniform_real_distribution<Double> dist;
	public:
		PhaseSpace() = delete; // must provide a number
		
		// initializes random number generator and arrays with the number of dimensions
		// for this particular phase space generator
		PhaseSpace(UInt32 _numDims) :
			gen(), dist(0.0, 1.0), numDims(_numDims),
			min(new Double[_numDims]), max(new Double[_numDims]), delta(new Double[_numDims])
		{}
		~PhaseSpace();

		inline const UInt32  getDim()    { return numDims; }
		inline const Double* getMins()   { return min;     }
		inline const Double* getMaxes()  { return max;     }
	    inline const Double* getDeltas() { return delta;   }

	protected:
		inline void fillDelta() {
			for (int i=0; i<numDims; i++)
				delta[i] = max[i] - min[i];
		}
	};


	// specifically considers the phase space specified by
	// cos(theta), y (the rapidity), and tau (determined by rho, technically,
	// which is actually what we generate a random number for)
	class PhaseSpace_TauYCosth : public PhaseSpace {
	private:
		// these are determined via the ECM
		Double rhoMin;
		Double rhoMax;
		Double deltaRho;
	public:
		PhaseSpace_TauYCosth();
		~PhaseSpace_TauYCosth() = default;
	};


	class PhaseSpace_Costh : public PhaseSpace {
	public:
		PhaseSpace_Costh();
		~PhaseSpace_Costh() = default;
	};
};


#endif // __PHASE_SPACE_HPP
