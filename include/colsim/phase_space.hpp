#ifndef __PHASE_SPACE_HPP
#define __PHASE_SPACE_HPP

#include "colsim/common.hpp"

#include <vector>
#include <initializer_list>

namespace colsim
{

	class PhaseSpace
	{
	protected:
		uint _num_dims;
		std::vector<double> _min, _max, _delta;

		// --- to be used for plotting ---
		// string names for phase space elements
		std::vector<std::string> _names;
		std::vector<std::string> _titles;
		std::vector<std::string> _xlabels;
		std::vector<std::string> _ylabels;
		
	public:
		PhaseSpace() = delete; // must provide a dimension
		
		// initializes random number generator and arrays with the number of dimensions
		// for this particular phase space generator
		PhaseSpace(uint num_dims,
				   std::initializer_list<std::string> names,
				   std::initializer_list<std::string> titles={},
				   std::initializer_list<std::string> xlabels={},
				   std::initializer_list<std::string> ylabels={}) :
			_num_dims{num_dims},
			_names{names}, _titles{titles}, _xlabels{xlabels}, _ylabels{ylabels}
		{}
		~PhaseSpace() = default;

		// some getters
		inline uint  dims() const { return _num_dims; }
		inline std::vector<double> const& mins()       const { return _min;   }
		inline std::vector<double> const& maxes()      const { return _max;   }
	    inline std::vector<double> const& deltas()     const { return _delta; }

		inline std::vector<std::string> const& names()   const { return _names; }
		inline std::vector<std::string> const& titles()  const { return _titles; }
		inline std::vector<std::string> const& xlabels() const { return _xlabels; }
		inline std::vector<std::string> const& ylabels() const { return _ylabels; }
		
		/** Fills @a vec with @a numDims elements corresponding
		 *  to randomly generated phase space points.
		 */
		void fill_phase_space(std::vector<double>& vec);

	protected:
		/** simple helper for filling the delta values
		 *  once the mins and maxes have been provided
		 */
		void fill_delta();


	public:
		/** Main routine that sets mins, maxes, deltas.
		 *  called by the constructor as well as when
		 *  any config file variables change
		 */
		void set_ranges() { }
	};


	// specifically considers the phase space specified by
	// cos(theta), y (the rapidity), and tau (determined by rho, technically,
	// which is actually what we generate a random number for)
	class PhaseSpace_TauYCosth : public PhaseSpace
	{
	private:
		// these are determined via the ECM
		double _rho_min;
		double _rho_max;
		double _drho;
	public:
		PhaseSpace_TauYCosth();
		~PhaseSpace_TauYCosth() = default;
		
		void set_ranges();
	};


	class PhaseSpace_EtEta : public PhaseSpace
	{
	private:
		
	public:
		PhaseSpace_EtEta();
		~PhaseSpace_EtEta() = default;
		
		void set_ranges();
	};
};


#endif // __PHASE_SPACE_HPP
