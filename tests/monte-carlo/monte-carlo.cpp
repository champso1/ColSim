#include "ColSim/ColSim.hpp"
#include <ColSim/Math.hpp>
using namespace ColSim;
using namespace std;

#include <iostream>
#include <fstream>

double testFunc1(double x[]) {
	return 2*x[0];
}

double testFunc2(double x[]) {
	return 4.0/(1+x[0]*x[0]);
}

double testFunc3(double x[]) {
	return exp(-x[0]*x[0] - x[1]*x[1] - x[2]*x[2]);
}

const static double ACTUAL1 = 1.0;
const static double ACTUAL2 = M_PI;
const static double ACTUAL3 = pow(M_PI, 1.5);

int main() {
	ColSimMain colsim;
	// colsim.init(ColSimMain::HARD_SCATTERING, "../../../res/config.in");
	colsim.stop();

	{
		MonteCarloParams params;
		params.func = testFunc1;
		double min[] = {0.0}, max[] = {1.0};
		params.min = min;
		params.max = max;
		params.numDims = 1.0;
		params.numEvals = 1000000;
		IntegrationResult res = MonteCarloIntegrate(params);

		double relError = abs(ACTUAL1-res.result)/ACTUAL1;

		LOGGER.logMessage("Calculated: %.7lf +- %.7lf", res.result, res.error);
		LOGGER.logMessage("Expected: %.1lf", ACTUAL1);
		LOGGER.logMessage("Relative Error: %.7lf", relError);

		if (relError >= 0.01) {
			LOGGER.logAbort("Monte Carlo algorithm has failed to converge for a sufficiently high number of iterations.");
		}
	}

	{
		MonteCarloParams params;
		params.func = testFunc2;
		double min[] = {0.0}, max[] = {1.0};
		params.min = min;
		params.max = max;
		params.numDims = 1.0;
		params.numEvals = 1000000;
		IntegrationResult res = MonteCarloIntegrate(params);

		double relError = abs(ACTUAL2-res.result)/ACTUAL2;

		LOGGER.logMessage("Calculated: %.7lf +- %.7lf", res.result, res.error);
		LOGGER.logMessage("Expected: %.7lf", ACTUAL2);
		LOGGER.logMessage("Relative Error: %.7lf", relError);

		if (relError >= 0.01) {
			LOGGER.logAbort("Monte Carlo algorithm has failed to converge for a sufficiently high number of iterations.");
		}
	}

	{
		MonteCarloParams params;
		params.func = testFunc3;
		double min[] = {-100.0, -100.0, -100.0}, max[] = {100.0, 100.0, 100.0};
		params.min = min;
		params.max = max;
		params.numDims = 3.0;
		params.numEvals = 1000000000;
		IntegrationResult res = MonteCarloIntegrate(params);

		double relError = abs(ACTUAL3-res.result)/ACTUAL2;

		LOGGER.logMessage("Calculated: %.7lf +- %.7lf", res.result, res.error);
		LOGGER.logMessage("Expected: %.7lf", ACTUAL3);
		LOGGER.logMessage("Relative Error: %.7lf", relError);

		if (relError >= 0.1) {
			LOGGER.logAbort("Monte Carlo algorithm has failed to converge for a sufficiently high number of iterations.");
		}
	}
	
   	return 0;
}
