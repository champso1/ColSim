#include "colsim/math.hpp"
#include "colsim/utils.hpp"

#include <cmath>

namespace colsim
{
	int rand_int(uint low, uint high)
	{
		static std::uniform_int_distribution<> dist(low, high);
		return dist(random_utils.mt);
	}

	double rand_double(double low, double high)
	{
		static std::uniform_real_distribution<double> dist(low, high);
		return dist(random_utils.mt);
	}

	double rand_double()
	{
		static std::uniform_real_distribution<double> dist(0.0, 1.0);
		return dist(random_utils.mt);
	}

}; // namespace ColSim

