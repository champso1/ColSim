#ifndef __PHASE_SPACE_HPP
#define __PHASE_SPACE_HPP

#include "ColSim/Types.hpp"
#include "ColSim/Constants.hpp"

#include <array>
#include <initializer_list>
#include <random>

namespace ColSim {

	class PhaseSpace {
	protected:
		USize numDims;
		std::vector<Double> min, max, delta;

		// --- to be used for plotting ---
		// string names for phase space elements
		std::vector<std::string> names;
		std::vector<std::string> titles;
		std::vector<std::string> xlabels;
		std::vector<std::string> ylabels;
		
	public:
		PhaseSpace() = delete; // must provide a dimension
		
		// initializes random number generator and arrays with the number of dimensions
		// for this particular phase space generator
		PhaseSpace(UInt32 _numDims,
				   std::initializer_list<std::string> _names,
				   std::initializer_list<std::string> _titles={},
				   std::initializer_list<std::string> _xlabels={},
				   std::initializer_list<std::string> _ylabels={}) :
			numDims(_numDims),
			names(_names), titles(_titles), xlabels(_xlabels), ylabels(_ylabels)
		{}
		~PhaseSpace() = default;

		// some getters
		inline UInt32  getDim() const { return numDims; }
		inline const std::vector<Double>& getMins()       const { return min;   }
		inline const std::vector<Double>& getMaxes()      const { return max;   }
	    inline const std::vector<Double>& getDeltas()     const { return delta; }

		inline const std::vector<std::string>& getNames()   const { return names; }
		inline const std::vector<std::string>& getTitles()  const { return titles; }
		inline const std::vector<std::string>& getXLabels() const { return xlabels; }
		inline const std::vector<std::string>& getYLabels() const { return ylabels; }
		
		/** Fills @a vec with @a numDims elements corresponding
		 *  to randomly generated phase space points.
		 */
		void fillPhaseSpace(std::vector<Double>& vec);

	protected:
		/** simple helper for filling the delta values
		 *  once the mins and maxes have been provided
		 */
		void fillDelta();


	public:
		/** Main routine that sets mins, maxes, deltas.
		 *  called by the constructor as well as when
		 *  any config file variables change
		 */
		void setRanges() { }
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
		
		void setRanges();
	};


	class PhaseSpace_EtEta : public PhaseSpace {
	private:
		
	public:
		PhaseSpace_EtEta();
		~PhaseSpace_EtEta() = default;
		
		void setRanges();
	};
};


#endif // __PHASE_SPACE_HPP
