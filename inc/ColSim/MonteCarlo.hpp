#ifndef __MONTE_CARLO_HPP
#define __MONTE_CARLO_HPP

#include "ColSim/Types.hpp"

#include <iostream>
#include <iomanip>
#include <limits>
#include <functional>


namespace ColSim {

#define DEFAULT_SEED 12345

	struct IntegrationResult {
		Double result;
		Double error, variance;
		Double maxWeight;
		USize numDims;
		Double* maxPoints;

		IntegrationResult(const USize _numDims)
			: result(0.0F), error(0.0F), variance(0.0F), maxWeight(0.0F),
			  numDims(_numDims), maxPoints(new Double[_numDims])
		{}
		~IntegrationResult() { delete[] maxPoints; }

		// for ease of printing to some output stream
		friend std::ostream& operator<<(std::ostream& s,
										const IntegrationResult& self)
		{
			s << std::setprecision(std::numeric_limits<double>::max_digits10)
			  << self.result << " +- " << self.error << "\n";
			return s;
		}
	};

	// follows the same form of the GSL integration routines
	// the only parameter is an array of the independent variables
	// unlike gsl, there is no second void* param,
	// since those we get from the Settings global
	typedef Double (*MonteCarloFunc)(Double*);


	struct MonteCarloParams {
		USize num_evals;
		USize num_dims;
		const Double* min;
		const Double* max;
		std::function<Double(Double*)> func;
	};



	extern IntegrationResult MonteCarloIntegrate(const MonteCarloParams& params);
}


#endif // __MONTE_CARLO_HPP
