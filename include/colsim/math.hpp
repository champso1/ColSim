#ifndef __MATH_HPP
#define __MATH_HPP

#include <limits>
#include <random>
#include <type_traits>

#include "colsim/utils.hpp"

namespace colsim
{
	template <typename T>
	concept MathFunction =
		std::is_invocable_v<T, double, void*> &&
		std::same_as<std::invoke_result_t<T, double, void*>, double>;

	struct RandomUtils final
	{
		std::random_device dev;
		std::mt19937_64 mt;

		RandomUtils() : dev{}, mt(dev()) {}
	};
    int rand_int(uint low, uint high);
	double rand_double(double low, double high);
	double rand_double();
	inline RandomUtils random_utils{};

	// WARNING: this requires an initial guess,
	// something which is very hard to accurately do,
	// so just do the bisection method
	template <MathFunction TMathFunction>
	double NewtonRaphson(TMathFunction func, double x0, uint numEvals, double precision, double dx, void* params)
	{
		double n = 0;
		double val = std::numeric_limits<double>::max();

		while ((std::abs(val) > precision) && (n < numEvals)) {
			double CD = (func(x0+dx/2.0, params) - func(x0-dx/2.0, params))/dx;
			double Dx = -func(x0,params)/CD;

			x0 += Dx;
			val = func(x0,params);
			n += 1;
		}

		if (n >= numEvals)
			log(LOG_ERROR, "NewtonRaphson():", "Failed to converge.");
			
	    return x0;
	}

	template <MathFunction TMathFunction>
	double Bisection(TMathFunction func, double min, double max, uint numEvals, double precision, void* params)
	{
		int count = 0;
		double range = max - min;
		double midpoint = std::numeric_limits<double>::max(); // arbitrary large number to obviously indicate error
		
		while ((range > precision) && (count < numEvals)) {
			midpoint = (min+max)/2.0;

			double funcMin = func(min, params);
			double funcMax = func(max, params);
			double funcMidpoint = func(midpoint, params);

			if ((funcMin * funcMidpoint) < 0)
				max = midpoint;
			else if ((funcMidpoint * funcMax) < 0)
				min = midpoint;
			else
				break;

			count += 1;
			range = max - min;
		}

		return midpoint;
	}
	

	// structure to hold information obtained as a result of
	// Monte Carlo integration
	struct IntegrationResult final
	{
		uint numDims;
		double result;
		double error, variance;

		// we are interested in the maximum values
		// achieved by the weight and each
		// of the phase space points
		double maxWeight;
		std::vector<double> maxPoints;


		// requires nothing special
		IntegrationResult(uint dims)
			: numDims{dims}, maxPoints(dims, 0.0)
		{}
		~IntegrationResult() = default;
	};

	// roughly follows the same form of the GSL integration routines
	// but without the void* since all such data
	// is in the settings singleton
	template <typename T>
	concept MonteCarloFunction =
		std::is_invocable_v<T, double[]> &&
		std::same_as<std::invoke_result_t<T, double[]>, double>;

	template <MonteCarloFunction TMonteCarloFunction>
	struct MonteCarloParams final
	{
		using func_type = TMonteCarloFunction;

		uint numEvals;
		uint numDims;
		double const* min;
		double const* max;
		func_type func;
	};
	
	template <MonteCarloFunction TMonteCarloFunction>
	IntegrationResult MonteCarloIntegrate(MonteCarloParams<TMonteCarloFunction> const& params)
	{
		double weightSum = 0.0F;
		double weightSquaredSum = 0.0F;
		IntegrationResult result(params.numDims);
		double points[params.numDims];
		
		// here we calculate the delta values
		// since we multiply by the result of the function
		// by definition of monte carlo integration
		double deltaX[params.numDims];
		// Double* delta_x = new Double[params.numDims];
		for (uint i=0; i<params.numDims; i++)
			deltaX[i] = params.max[i] - params.min[i];

		
		for (uint i=0; i<params.numEvals; i++) {
			// grab the random value for each dim
			for (uint j=0; j<params.numDims; j++) 
			    points[j] = params.min[j] + rand_double() * deltaX[j];
			
			// calculate that sheisse
			double weight = (params.func)(points);
	
			// multiply by the deltas of the independent variables
			for (uint j=0; j<params.numDims; j++)
			    weight *= deltaX[j];

			// add to total weight
			weightSum += weight;
			weightSquaredSum += pow(weight, 2);

			// set max stuff
			if (weight > result.maxWeight) {
				result.maxWeight = weight;
				for (uint j=0; j<params.numDims; j++)
					result.maxPoints[j] = points[j];
			}
		}
		
		double numEvals_d = static_cast<double>(params.numEvals);
		
		result.result = weightSum/numEvals_d;
		result.variance = weightSquaredSum/numEvals_d - std::pow(weightSum/numEvals_d,2);
		result.error = std::sqrt(result.variance/numEvals_d);

		return result;
	}
}


#endif // __MATH_HPP
