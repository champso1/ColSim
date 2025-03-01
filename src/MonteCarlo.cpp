#include "ColSim/MonteCarlo.hpp"

#include <random>
#include <cmath>
#include <functional>

namespace ColSim {

    IntegrationResult MonteCarloIntegrate(const MonteCarloParams& params) {
		std::default_random_engine engine;
		engine.seed(DEFAULT_SEED);
		std::uniform_real_distribution<Double> distribution(0.0, 1.0F);
		auto rand = std::bind(distribution, engine);
		
		Double weightSum = 0.0F;
		Double weightSquaredSum = 0.0F;

		IntegrationResult result(params.num_dims);
		
		// array for storing the random values
		// to be passed to the function
		Double* points = new Double[params.num_dims];

		// here we calculate the delta values
		// since we multiply by the result of the function
		// by definition of monte carlo integration
		Double* delta_x = new Double[params.num_dims];
		for (int i=0; i<params.num_dims; i++) 
			delta_x[i] = params.max[i] - params.min[i];

		
		for (int i=0; i<params.num_evals; i++) {
			// grab the random value for each dim
			for (int j=0; j<params.num_dims; j++) 
			    points[j] = params.min[j] + rand() * delta_x[j];
			
			// calculate that sheisse
			Double weight = (params.func)(points);
	
			// multiply by the deltas of the independent variables
			for (int j=0; j<params.num_dims; j++)
			    weight *= delta_x[j];

			// add to total weight
			weightSum += weight;
			weightSquaredSum += pow(weight, 2);

			// set max stuff
			if (weight > result.maxWeight) {
				result.maxWeight = weight;
				for (int j=0; j<params.num_dims; j++)
					result.maxPoints[j] = points[j];
			}
		}
		delete[] points;
		delete[] delta_x;

		Double num_evals_d = static_cast<Double>(params.num_evals);
		result.result = weightSum / num_evals_d;

		result.variance = (weightSquaredSum / num_evals_d - pow(weightSum / num_evals_d, 2));
		result.error = sqrt(result.variance / num_evals_d);

		return result;
	}


}; // namespace ColSim

