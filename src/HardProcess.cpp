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
		std::vector<Double> points;

		// here we get delta values which we need for Monte Carlo
	    const std::vector<Double>& deltas = phaseSpace->getDeltas();
		
		for (UInt i=0; i<SETTINGS.numXSIterations; i++) {
			phaseSpace->fillPhaseSpace(points);
			
			// calculate that sheisse
			Result dsigmaRes = dSigma(points);

			// if it is invalid, we must redo
			if (dsigmaRes == Result::invalidResult()) {
				i--;
				continue;
			}
			
			Double weight = dsigmaRes.weight;
    
			// multiply by the deltas of the independent variables
			for (UInt j=0; j<numDims; j++)
			    weight *= deltas[j];

			// add to total weight
			weightSum += weight;
			weightSquaredSum += pow(weight, 2);

			// set maxes
			if (weight > res.maxWeight) {
				res.maxWeight = weight;
				for (UInt j=0; j<numDims; j++)
					res.maxPoints[j] = points[j];
			}
		}

		// doing divisions, so we want a Double
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


	HardProcess::Result PP2Zg2ll::dSigma(const std::vector<Double>& phaseSpacePoints) {
		const Double S = SETTINGS.S; // ECM^2
		const double E_TR_2 = SETTINGS.transEnergy_2;

		// independent variables
		const Double cosTheta = phaseSpacePoints[0];
		const Double rho      = phaseSpacePoints[1];
		const double rand_y   = phaseSpacePoints[2];

		// other variables
	    const Double jacobian = (E_TR_2) / (std::cos(rho)*std::cos(rho) * S);

		const Double s_hat = E_TR_2*std::tan(rho) + E_TR_2;

		const Double yMax   = -0.5*std::log(s_hat/S);
		const Double deltaY = 2.0*yMax;
		const Double y      = (2.0*rand_y - 1.0)*yMax;

		const Double x1 = std::sqrt(s_hat/S)*std::exp(y);
		const Double x2 = std::sqrt(s_hat/S)*std::exp(-y);

		if ((x1 > 1.0 || x1 < 0.0) || (x2 > 1.0 || x2 <0.0)) {
			LOGGER.logWarning("Encountered invalid x1 or x2: (%.4lf,%.4lf). Retrying...", x1, x2);
			return HardProcess::Result::invalidResult();
		}

		Double weight = computeWeight(s_hat, x1, x2, cosTheta);
		// go ahead and scale it this way here
		weight *= (jacobian * deltaY);
		weight /= (x1 * x2);

		// pass through the COM energy and the momentum fraction (x1)
		return HardProcess::Result(weight, {std::sqrt(s_hat), x1});
	}


	void PP2Zg2ll::generateParticles(std::vector<Particle>& particles) {
		const Double S = SETTINGS.S;
		const Double ECM = SETTINGS.ECM;
		const double E_TR_2 = SETTINGS.transEnergy_2;

		std::vector<Double> phaseSpacePoints;
		phaseSpace->fillPhaseSpace(phaseSpacePoints);

		const Double cosTheta = phaseSpacePoints[0];
		const Double rho      = phaseSpacePoints[1];
		const double rand_y   = phaseSpacePoints[2];

		// other variables
		const Double s_hat = E_TR_2*std::tan(rho) + E_TR_2;
		const Double Q = std::sqrt(s_hat);

		const Double yMax   = -0.5*std::log(s_hat/S);
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












	PP2Jets::PP2Jets() {
		phaseSpace = std::unique_ptr<PhaseSpace>(new PhaseSpace_EtEta());
	}


	Double PP2Jets::qqprime2qqprime(Double eta3, Double eta4, Double E_t, Double alpha_s) {
		Double eta_star = calcEtaStar(eta3, eta4);
		Double s = calcSHat(eta_star, E_t);
		Double t = calcTHat(eta_star, E_t);
		Double u = calcUHat(eta_star, E_t);

		return (4*M_PI*alpha_s*alpha_s)/(9*s*s) * ((s*s +  u*u)/(t*t));
	}
	Double PP2Jets::qqprimebar2qqprimebar(Double eta3, Double eta4, Double E_t, Double alpha_s) {
		return qqprime2qqprime(eta3, eta4, E_t, alpha_s);
	}

	// scattering of identical quarks
	Double PP2Jets::qq2qq(Double eta3, Double eta4, Double E_t, Double alpha_s) {
		Double eta_star = calcEtaStar(eta3, eta4);
		Double s = calcSHat(eta_star, E_t);
		Double t = calcTHat(eta_star, E_t);
		Double u = calcUHat(eta_star, E_t);

		return (4.0*M_PI*alpha_s*alpha_s)/(9.0*s*s) * ((s*s + u*u)/(t*t) + (t*t + u*u)/(s*s) - (2.0/3.0)*((s*s)/(u*t)));
	}

	Double PP2Jets::qqbar2qprimeqprimebar(Double eta3, Double eta4, Double E_t, Double alpha_s) {
		Double eta_star = calcEtaStar(eta3, eta4);
		Double s = calcSHat(eta_star, E_t);
		Double t = calcTHat(eta_star, E_t);
		Double u = calcUHat(eta_star, E_t);


		return (4.0*M_PI*alpha_s*alpha_s)/(9.0*s*s) * ((t*t + u*u)/(s*s));
	}


	Double PP2Jets::gg2gg(Double eta3, Double eta4, Double E_t, Double alpha_s) {
		Double eta_star = calcEtaStar(eta3, eta4);
		Double s = calcSHat(eta_star, E_t);
		Double t = calcTHat(eta_star, E_t);
		Double u = calcUHat(eta_star, E_t);

		return (9.0*M_PI*alpha_s*alpha_s) / (2.0 * s*s) * (3 - (t*u)/(s*s) - (s*u)/(t*t) - (s*t)/(u*u));
    
	}
	Double PP2Jets::qqbar2gg(Double eta3, Double eta4, Double E_t, Double alpha_s) {
		Double eta_star = calcEtaStar(eta3, eta4);
		Double s = calcSHat(eta_star, E_t);
		Double t = calcTHat(eta_star, E_t);
		Double u = calcUHat(eta_star, E_t);

		return (32.0*M_PI*alpha_s*alpha_s)/(27.0*s*s) * (u/t + t/u - (9.0/4.0)*((t*t + u*u)/(s*s)));
	}
	Double PP2Jets::gg2qqbar(Double eta3, Double eta4, Double E_t, Double alpha_s) {
		Double eta_star = calcEtaStar(eta3, eta4);
		Double s = calcSHat(eta_star, E_t);
		Double t = calcTHat(eta_star, E_t);
		Double u = calcUHat(eta_star, E_t);

		return (M_PI * alpha_s*alpha_s)/(6*s*s) * (u/t + t/u - (9.0/4.0)*((t*t + u*u)/(s*s)));
	}

	Double PP2Jets::qg2qg(Double eta3, Double eta4, Double E_t, Double alpha_s) {
		Double eta_star = calcEtaStar(eta3, eta4);
		Double s = calcSHat(eta_star, E_t);
		Double t = calcTHat(eta_star, E_t);
		Double u = calcUHat(eta_star, E_t);

		return (4.0*M_PI*alpha_s*alpha_s)/(9*s*s) * (-u/s - s/u + (9.0/4.0)*((s*s + u*u)/(t*t)));
	}


	HardProcess::Result PP2Jets::dSigma(const std::vector<Double>& phaseSpacePoints) {
		const Double Et = phaseSpacePoints[0];
		const Double eta3 = phaseSpacePoints[1];

		// retrieve/calculate quantities required to grab PDF data
		const Double S = SETTINGS.S;
		LHAPDF::PDF* pdf = SETTINGS.pdf;
		Double alphaS = pdf->alphasQ2(S);
		Double x1 = calcX1(eta3, _eta4, Et);
		Double x2 = calcX2(eta3, _eta4, Et);
		const std::vector<UInt32> pids{1, 2, 3, 4, 5, 21}; // vector of PDG IDs for all quarks (but top) and gluon

		// must add all processes
		Double total = 0.0;

		// q + q' -> q + q'
		Double total_qqprime2qqprime = 0.0;
		for (UInt i=0; i<pids.size(); i++) {
			for (int j=i+1; i<pids.size(); i++) {
				total_qqprime2qqprime += pdf->xfxQ2(pids[i], x1, S)*pdf->xfxQ2(pids[j], x2, S)
					* (1.0/M_PI)*qqprime2qqprime(eta3, _eta4, Et, alphaS);
			}
		}
		total += total_qqprime2qqprime;

		// q + q -> q + q
		Double total_qq2qq = 0.0;
		for (UInt i=0; i<pids.size(); i++) {
			total_qqprime2qqprime += pdf->xfxQ2(pids[i], x1, S)*pdf->xfxQ2(pids[i], x2, S)
				* (1.0/M_PI)*qq2qq(eta3, _eta4, Et, alphaS);
		}
		total += total_qq2qq;

		// q + qb -> q' + qb'
		Double total_qqbar2qprimeqprimebar = 0.0;
		for (UInt i=0; i<pids.size(); i++) {
			total_qqprime2qqprime += pdf->xfxQ2(pids[i], x1, S)*pdf->xfxQ2(pids[i], x2, S)
				* (1.0/M_PI)*qqbar2qprimeqprimebar(eta3, _eta4, Et, alphaS);
		}
		total += total_qqbar2qprimeqprimebar;

		
		// g + g -> g + g
	    total += pdf->xfxQ2(21, x1, S) * pdf->xfxQ2(21, x2, S) * gg2gg(eta3, _eta4, Et, alphaS);

		// q + qb -> g + g
		Double total_qqbar2gg = 0.0;
		for (UInt i=0; i<pids.size(); i++) {
			total_qqprime2qqprime += pdf->xfxQ2(pids[i], x1, S)*pdf->xfxQ2(pids[i], x2, S)
				* (1.0/M_PI)*qqbar2gg(eta3, _eta4, Et, alphaS);
		}
		total += total_qqbar2gg;

		// g + g -> q + qb
		Double total_gg2qqbar = 0.0;
		for (UInt i=0; i<pids.size(); i++) {
			total_qqprime2qqprime += pdf->xfxQ2(pids[i], x1, S)*pdf->xfxQ2(pids[i], x2, S)
				* (1.0/M_PI)*gg2qqbar(eta3, _eta4, Et, alphaS);
		}
		total += total_gg2qqbar;

		// q + g -> q + g
		Double total_qg2qg = 0.0;
		for (UInt i=0; i<pids.size(); i++) {
			total_qqprime2qqprime += pdf->xfxQ2(pids[i], x1, S)*pdf->xfxQ2(pids[i], x2, S)
				* (1.0/M_PI)*qg2qg(eta3, _eta4, Et, alphaS);
		}
		total += 2*total_qg2qg;

		total *= 2.0*Et;
		return Result(total);
	}

	void PP2Jets::generateParticles(std::vector<Particle>& momenta) {
		// does nothing yet
		(void)momenta;
		return;
	}

}; // namespace ColSim
