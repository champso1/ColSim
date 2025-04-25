#ifndef __MATH_HPP
#define __MATH_HPP

#include "ColSim/Types.hpp"

#include <iostream>
#include <iomanip>
#include <limits>
#include <functional>
#include <random>

namespace ColSim {

	// random stuff
#define DEFAULT_RANDOM_SEED 12345

	extern std::default_random_engine engine;
	extern std::uniform_real_distribution<double> distribution;
	extern std::function<Double()> randDouble;

	typedef std::function<Double(Double,void*)> MathFunction;

	// WARNING: this requires an initial guess,
	// something which is very hard to accurately do,
	// so just do the bisection method
	extern Double NewtonRaphson(MathFunction func, Double x0, UInt32 numEvals, Double precision, Double dx, void* params);
	extern Double Bisection(MathFunction func, Double min, Double max, UInt32 numEvals, Double precision, void* params);
	

	// structure to hold information obtained as a result of
	// Monte Carlo integration
	struct IntegrationResult {
		USize numDims;
		
		Double result;
		Double error, variance;

		// we are interested in the maximum values
		// achieved by the weight and each
		// of the phase space points
		Double maxWeight;
		std::vector<Double> maxPoints;


		// requires nothing special
		IntegrationResult(const USize _numDims)
			: numDims(_numDims), maxPoints(_numDims, 0.0)
		{}
		~IntegrationResult() = default;
	};

	// roughly follows the same form of the GSL integration routines
	// but without the void* since all such data
	// is in the settings singleton
	typedef std::function<Double(Double[])> MonteCarloFunc;

	struct MonteCarloParams {
		USize numEvals;
		USize numDims;
		const Double* min;
		const Double* max;
		MonteCarloFunc func;
	};

	/** Performs a monte carlo integration with the function given in @a params.
	 *
	 *  @returns An @a IntegrationResult that contains useful information
	 *   (including the result!) from the integration.
	 */
	extern IntegrationResult MonteCarloIntegrate(const MonteCarloParams& params);
}


#endif // __MATH_HPP
