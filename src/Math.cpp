#include "ColSim/Math.hpp"
#include "ColSim/Logger.hpp"

#include <cmath>
#include <cfloat>

namespace ColSim {

	std::default_random_engine engine;
	std::uniform_real_distribution<double> distribution(0, 1.0);
	std::function<Double()> randDouble;

	Double NewtonRaphson(MathFunction func, Double x0, UInt32 numEvals, Double precision, Double dx, void* params) {
		Double n = 0;
		Double val = DBL_MAX;

		LOGGER.logMessage("Converging via Newton-Raphson...");
		
		while ((std::abs(val) > precision) && (n < numEvals)) {
			Double CD = (func(x0+dx/2.0, params) - func(x0-dx/2.0, params))/dx;
			Double Dx = -func(x0,params)/CD;

			x0 += Dx;
			val = func(x0,params);
			n += 1;

			LOGGER.logMessage("\tx0 -> %.9lf", x0);
		}

		if (n >= numEvals) 
			throw std::runtime_error("Failed to converge Newton Raphson root-finding.\n");

	    return x0;
	}


	Double Bisection(MathFunction func, Double min, Double max, UInt32 numEvals, Double precision, void* params) {
		UInt32 count = 0;
		Double range = max - min;
		Double midpoint = DBL_MAX; // arbitrary large number to obviously indicate error
		
		LOGGER.logMessage("Converging via bisection...");
		
		while ((range > precision) && (count < numEvals)) {
			midpoint = (min+max)/2.0;

			Double funcMin = func(min, params);
			Double funcMax = func(max, params);
			Double funcMidpoint = func(midpoint, params);

			if ((funcMin * funcMidpoint) < 0)
				max = midpoint;
			else if ((funcMidpoint * funcMax) < 0)
				min = midpoint;
			else
				break;

			count += 1;
			range = max - min;

		    LOGGER.logMessage("\t%u [%.9lf, %.9lf]", count, min, max);
		}

		return midpoint;
	}
	

    IntegrationResult MonteCarloIntegrate(const MonteCarloParams& params) {
		Double weightSum = 0.0F;
		Double weightSquaredSum = 0.0F;

		IntegrationResult result(params.numDims);
		
		// array for storing the random values
		// to be passed to the function
		Double points[params.numDims];
		// Double* points = new Double[params.numDims];

		// here we calculate the delta values
		// since we multiply by the result of the function
		// by definition of monte carlo integration
		Double deltaX[params.numDims];
		// Double* delta_x = new Double[params.numDims];
		for (int i=0; i<params.numDims; i++) 
			deltaX[i] = params.max[i] - params.min[i];

		
		for (int i=0; i<params.numEvals; i++) {
			// grab the random value for each dim
			for (int j=0; j<params.numDims; j++) 
			    points[j] = params.min[j] + randDouble() * deltaX[j];
			
			// calculate that sheisse
			Double weight = (params.func)(points);
	
			// multiply by the deltas of the independent variables
			for (int j=0; j<params.numDims; j++)
			    weight *= deltaX[j];

			// add to total weight
			weightSum += weight;
			weightSquaredSum += pow(weight, 2);

			// set max stuff
			if (weight > result.maxWeight) {
				result.maxWeight = weight;
				for (int j=0; j<params.numDims; j++)
					result.maxPoints[j] = points[j];
			}
		}
		
		Double numEvals_d = static_cast<Double>(params.numEvals);
		
		result.result = weightSum/numEvals_d;
		result.variance = weightSquaredSum/numEvals_d
			              - std::pow(weightSum/numEvals_d,2);
		result.error = std::sqrt(result.variance/numEvals_d);

		return result;
	}


}; // namespace ColSim

