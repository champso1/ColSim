#ifndef __PHASE_SPACE_HPP
#define __PHASE_SPACE_HPP

#include "ColSim/Types.hpp"
#include "ColSim/Constants.hpp"

#include <array>
#include <random>

namespace ColSim {

	class PhaseSpace {
	protected:
		USize numDims;
		std::vector<Double> min, max, delta;
	    
	public:
		PhaseSpace() = delete; // must provide a dimension
		
		// initializes random number generator and arrays with the number of dimensions
		// for this particular phase space generator
		PhaseSpace(UInt32 _numDims) :
			numDims(_numDims),
			min(_numDims), max(_numDims), delta(_numDims)
		{}
		~PhaseSpace() = default;

		// some getters
		inline const UInt32  getDim() const { return numDims; }
		inline const std::vector<Double>& getMins()   const { return min;   }
		inline const std::vector<Double>& getMaxes()  const { return max;   }
	    inline const std::vector<Double>& getDeltas() const { return delta; }

		/** Fills @a vec with @a numDims elements corresponding
		 *  to randomly generated phase space points.
		 */
		void fillPhaseSpace(std::vector<Double>& vec);
		

	protected:
		// simple helper for filling the delta values
		// once the mins and maxes have been provided
		void fillDelta();
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
};


#endif // __PHASE_SPACE_HPP
