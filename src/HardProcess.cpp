#include "ColSim/HardProcess.hpp"

#include "ColSim/Logger.hpp"
#include "ColSim/Settings.hpp"
#include "ColSim/Constants.hpp"
#include "ColSim/PhaseSpace.hpp"

#include "LHAPDF/LHAPDF.h"

#include <cmath>
#include <memory>

namespace ColSim {


	HardProcessResult HardProcess::calculate() {
		USize numDims = phaseSpace->getDim();
		
		Double weightSum = 0.0F;
		Double weightSquaredSum = 0.0F;

		HardProcessResult res(numDims);
		
	    // setup vector for phase space points
		std::vector<Double> points(numDims);

		// here we calculate the delta values
		// since we multiply by the result of the function
		// by definition of monte carlo integration
	    const std::vector<Double>& deltas = phaseSpace->getDeltas();
		
		for (int i=0; i<SETTINGS.numXSIterations; i++) {
			phaseSpace->fillPhaseSpace(points);
			
			// calculate that sheisse
			Double weight = dSigma(points);
	
			// multiply by the deltas of the independent variables
			for (int j=0; j<numDims; j++)
			    weight *= deltas[j];

			// add to total weight
			weightSum += weight;
			weightSquaredSum += pow(weight, 2);

			// set maxes
			if (weight > res.maxWeight) {
				res.maxWeight = weight;
				for (int j=0; j<numDims; j++)
					res.maxPoints[j] = points[j];
			}
		}
		
		Double numEvals_d = static_cast<Double>(SETTINGS.numXSIterations);
		
		res.result = weightSum/numEvals_d;
		Double variance = weightSquaredSum/numEvals_d
			              - std::pow(weightSum/numEvals_d,2);
		res.error = std::sqrt(variance/numEvals_d);

		// scale to picobarns
		res.result *= MAGIC_FACTOR;
		res.error *= MAGIC_FACTOR;
		
		return res;
	}



	

	PP2Zg2ll::PP2Zg2ll() {
		phaseSpace = std::unique_ptr<PhaseSpace>(new PhaseSpace_TauYCosth());
	}
	
	Double PP2Zg2ll::Kappa() const {
		return std::sqrt(2.0)*FERMI_CONSTANT*Z_MASS_2 / (4.0*PI*ALPHA);
	}

	Double PP2Zg2ll::Chi1(const double s_hat) const {
		return Kappa()*s_hat*(s_hat-Z_MASS_2) / (std::pow(s_hat-Z_MASS,2) + Z_WIDTH_2*Z_MASS_2);
	}
	
	Double PP2Zg2ll::Chi2(const Double s_hat) const {
		return std::pow(Kappa(),2) * pow(s_hat,2) / (std::pow(s_hat-Z_MASS_2,2) + Z_WIDTH_2*Z_MASS_2);
	}
	
	Double PP2Zg2ll::A0(const UInt8 quarkType, const Double s_hat) const {
		Double CAe = -0.5, CVe = -0.5 + 2*WEINBERG_ANGLE;
		Double CVf, CAf, Qf;
		if (quarkType == 0) { // up-type
			CVf = 0.5 - (4.0/3.0)*WEINBERG_ANGLE;
			CAf = 0.5;
			Qf = 2.0/3.0;
		} else { // down-type
			CVf = -0.5 + (2.0/3.0)*WEINBERG_ANGLE;
			CAf = -0.5;
			Qf = -1.0/3.0;
		}

	    return Qf*Qf - 2.0*Qf*CVe*CVf*Chi1(s_hat) + (CAe*CAe + CVe*CVe)*(CAf*CAf + CVf*CVf) * Chi2(s_hat);
	}
	
	Double PP2Zg2ll::A1(const UInt8 quarkType, const Double s_hat) const {
		Double A_mu = -0.5, V_mu = -0.5 + 2*WEINBERG_ANGLE;
		Double V_quark, A_quark, Q_quark;
		if (quarkType == 0) { // up-type
			V_quark = 0.5 -(4.0/3.0)*WEINBERG_ANGLE;
			A_quark = 0.5;
			Q_quark = 2.0/3.0;
		} else {
			V_quark = -0.5 + (2.0/3.0)*WEINBERG_ANGLE;
			A_quark = -0.5;
			Q_quark = -1.0/3.0;
		}

	    return -4.0*Q_quark*A_mu*A_quark*Chi1(s_hat) + 8.0*A_mu*V_mu*A_quark*V_quark*Chi2(s_hat);
	}


	Double PP2Zg2ll::dSigmaHat(const Double cosTheta, const UInt8 quarkType, const Double s_hat) {
		return (2.0*PI*std::pow(ALPHA,2) / (4.0*3.0*s_hat)) * (A0(quarkType,s_hat)*(1.0 + std::pow(cosTheta,2)) + A1(quarkType,s_hat)*cosTheta);
	}

	Double PP2Zg2ll::computeWeight(const Double s_hat, const Double x1, const Double x2, const Double cosTheta) {
		Double weight = 0.0F;
		LHAPDF::PDF* pdf = SETTINGS.pdf;
		// up-type quarks
		weight += dSigmaHat(cosTheta , 0, s_hat) * ((pdf->xfxQ2(2 , x1, s_hat) * pdf->xfxQ2(-2, x2, s_hat)) + (pdf->xfxQ2(4 , x1, s_hat) * pdf->xfxQ2(-4, x2, s_hat)));	
		weight += dSigmaHat(-cosTheta, 0, s_hat) * ((pdf->xfxQ2(-2, x1, s_hat) * pdf->xfxQ2(2 , x2, s_hat)) + (pdf->xfxQ2(-4, x1, s_hat) * pdf->xfxQ2(4 , x2, s_hat)));
		// down-type quarks
		weight += dSigmaHat(cosTheta , 1, s_hat) * ((pdf->xfxQ2(1 , x1, s_hat) * pdf->xfxQ2(-1, x2, s_hat)) + (pdf->xfxQ2(3 , x1, s_hat) * pdf->xfxQ2(-3, x2, s_hat)));
		weight += dSigmaHat(-cosTheta, 1, s_hat) * ((pdf->xfxQ2(-1, x1, s_hat) * pdf->xfxQ2(1 , x2, s_hat)) + (pdf->xfxQ2(-3, x1, s_hat) * pdf->xfxQ2(3 , x2, s_hat)));

		return weight;
	}


	Double PP2Zg2ll::dSigma(const std::vector<Double>& phaseSpacePoints) {
		LHAPDF::PDF* pdf = SETTINGS.pdf;
		const Double ECM = SETTINGS.ECM;
		const Double S = SETTINGS.S; // ECM^2

		// independent variables
		const Double cosTheta = phaseSpacePoints[0];
		const Double rho      = phaseSpacePoints[1];
		const double rand_y   = phaseSpacePoints[2];

		// other variables
	    const Double jacobian = (MASS_TR*WIDTH_TR) / (std::pow(std::cos(rho),2) * S);

		const Double s_hat = MASS_TR*WIDTH_TR*std::tan(rho) + std::pow(MASS_TR,2);

		const Double yMax   = -0.5*std::log(s_hat/S);
		const Double deltaY = 2.0*yMax;
		const Double y      = (2.0*rand_y - 1.0)*yMax;

		const Double x1 = std::sqrt(s_hat/S)*std::exp(y);
		const Double x2 = std::sqrt(s_hat/S)*std::exp(-y);

		Double weight = computeWeight(s_hat, x1, x2, cosTheta);
		// go ahead and scale it this way here
		weight *= (jacobian * deltaY);
		weight /= (x1 * x2);
		
		return weight;
	}


	void PP2Zg2ll::generateParticles(std::vector<Particle>& particles) {
		Double S = SETTINGS.S;
		Double ECM = SETTINGS.ECM;

		std::vector<Double> phaseSpacePoints;
		phaseSpace->fillPhaseSpace(phaseSpacePoints);

		const Double cosTheta = phaseSpacePoints[0];
		const Double rho      = phaseSpacePoints[1];
		const double rand_y   = phaseSpacePoints[2];

		// other variables
	    const Double jacobian = (MASS_TR*WIDTH_TR) / (std::pow(std::cos(rho),2) * S);

		const Double s_hat = MASS_TR*WIDTH_TR*std::tan(rho) + std::pow(MASS_TR,2);
		const Double Q = std::sqrt(s_hat);

		const Double yMax   = -0.5*std::log(s_hat/S);
		const Double deltaY = 2.0*yMax;
		const Double y      = (2.0*rand_y - 1.0)*yMax;

		const Double x1 = std::sqrt(s_hat/S)*std::exp(y);
		const Double x2 = std::sqrt(s_hat/S)*std::exp(-y);

		// boost parameter
		Double beta = (x2-x1)/(x2+x1);
		// more phase space points
		Double phi = randDouble()*2.0*M_PI;
		Double sinPhi = std::sin(phi);
		Double cosPhi = std::cos(phi);
		Double sinTheta = std::sqrt(1.0 - cosTheta*cosTheta);

		particles.emplace_back(FourVector(0.5*x1*ECM, 0.0, 0.0, 0.5*x1*ECM),
					1, "q1");
		particles.emplace_back(FourVector(0.5*x2*ECM, 0.0, 0.0, -0.5*x2*ECM),
					1, "q2");
		particles.emplace_back(FourVector(0.5*Q, 0.5*Q*sinTheta*cosPhi, 0.5*Q*sinTheta*sinPhi, 0.5*Q*cosTheta).zBoost(beta),
					13, "l1");
		particles.emplace_back(FourVector(0.5*Q, -0.5*Q*sinTheta*cosPhi, -0.5*Q*sinTheta*sinPhi, -0.5*Q*cosTheta).zBoost(beta),
					-13, "l2");
	}

}; // namespace ColSim
