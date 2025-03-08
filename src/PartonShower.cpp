#include "ColSim/PartonShower.hpp"
#include "ColSim/Constants.hpp"
#include "ColSim/Math.hpp"
#include "ColSim/Logger.hpp"
#include <fstream>

namespace ColSim {

	
	Double PartonShower::tGamma(Double z, Double alphaSOver) {
		return -2.0*alphaSOver*CF*std::log(1.0-z);
	}

	Double PartonShower::tGammaInverse(Double r, Double alphaSOver) {
		return 1.0 - std::exp(-0.5*r / (CF*alphaSOver));
	}

	Double PartonShower::zpOver(Double t, Double Qcut) {
		return 1.0-std::sqrt(Qcut*Qcut/t);
	}

	Double PartonShower::zmOver(Double t, Double Qcut) {
		return std::sqrt(Qcut*Qcut/t);
	}

	Double PartonShower::Pqq(Double z) {
		return CF*(1.0 + z*z) / (1.0-z);
	}
	
	Double PartonShower::PqqOver(Double z) {
		return 2.0*CF / (1.0-z);
	}


	Double PartonShower::emissionScaleFunc(Double log_T_Q2, void* _params) {
		Params* params = reinterpret_cast<Params*>(_params);
		Double Q = params->Q;
		Double Qcut = params->Qcut;
		Double r = params->r;
		Double alphaSOver = params->alphaSOver;

		Double t = Q*Q * std::exp(log_T_Q2);
		Double rho = tGamma(zpOver(t,Qcut), alphaSOver) - tGamma(zmOver(t,Qcut), alphaSOver);

		return log_T_Q2 - (1.0/rho)*std::log(r);
	}

    Double PartonShower::getTEmission(Double Q, Double Qcut, Double r, Double tfac, Double alphaSOver) {
		MathFunction emissionFunc = [this](Double x, void* y) -> Double {
			return this->emissionScaleFunc(x, y);
		};
		Params params{Q, Qcut, r, alphaSOver};

	    Double min = std::log(tfac * Qcut*Qcut / (Q*Q));
		Double max = 0.0f;
		Double tolerance = 1e-3F;

	    Double root = Bisection(emissionFunc, min, max, 1000, tolerance, reinterpret_cast<void*>(&params));

		Double tSol = Q*Q * std::exp(root);
		if (emissionScaleFunc(root, reinterpret_cast<void*>(&params)) > tolerance)
			return 0.0;
		return tSol;
	}



	Double PartonShower::getZEmission(Double t, Double Qcut, Double r, Double alphaSOver) {
	    return tGammaInverse( tGamma(zmOver(t, Qcut), alphaSOver) +
						  r * (tGamma(zpOver(t, Qcut), alphaSOver) - tGammaInverse(zmOver(t, Qcut), alphaSOver)), alphaSOver);
	}
	
	Double PartonShower::getPT_2(Double t, Double z) {
		return z*z * std::pow(1.0-z,2)*t;
	}
	
	Double PartonShower::getM_2(Double t, Double z) {
		return z*(1.0-z)*t;
	}

	PartonShower::Emission PartonShower::generateEmission(Double Q, Double Qcut, Double tfac, Double alphaSOver) {
		Bool generated = true;
		Double r1 = randDouble();
		Double r2 = randDouble();
		Double r3 = randDouble();
		Double r4 = randDouble();

		Double t = getTEmission(Q, Qcut, r1, tfac, alphaSOver);
		if (t == 0.0) 
			return Emission{t, 1.0, 0.0, 0.0, generated, false};
		

		Double z = getZEmission(t, Qcut, r2, alphaSOver);
		Double pT_2 = getPT_2(t, z);
		Double m_2 = getM_2(t, z);

		LOGGER.logMessage("Candidate emission scale: sqrt(t) = %.9lf", std::sqrt(t));
		LOGGER.logMessage("Candidate momentum fraction: z = %.9lf", z);
		LOGGER.logMessage("Candidate transverse momentum squared: pt^2 = %.9lf", pT_2);
		
		if (pT_2 < 0.0) {
			generated = false;
			LOGGER.logMessage("Emission rejected due to negative pT^2");
		}

		Double splittingFnOverRatio = Pqq(z) / PqqOver(z);
		if (splittingFnOverRatio < r3) {
		    LOGGER.logMessage("Emission rejected due to splitting function overestimate: p = %.9lf R = %.9lf", splittingFnOverRatio, r3);
			generated = false;
		} else {
			LOGGER.logMessage("Emission NOT rejected due to splitting function overestimate: p = %.9lf R = %.9lf", splittingFnOverRatio, r3);
		}

		Double alphaSOverRatio = alphaS.alphaSActual() / alphaS.alphaSOver();
		if (alphaSOverRatio < r4) {
			generated = false;
			LOGGER.logMessage("Emission rejected due to alpha_s overestimate: a_s = %.9lf R = %.9lf", alphaSOverRatio, r4);
		} else {
			LOGGER.logMessage("Emission NOT rejected due to alpha_s overestimate: a_s = %.9lf R = %.9lf", alphaSOverRatio, r4);
		}

		if (!generated) {
			z = 1.0;
			pT_2 = 0.0;
			m_2 = 0.0;
		}

		return Emission{t, z, pT_2, m_2, generated, true};
	}


	std::vector<PartonShower::Emission> PartonShower::Evolve(Double Q, Double Qmin, Double alphaSOver) {
		Double tMin = Qmin*Qmin;
		UInt32 numEmissions = 0;

		std::vector<Emission> emissions;

		Double fac_t = 3.9999;
		Double fac_cutoff = 4.0;

		Double t = Q*Q;
		Double z = 1.0;

		LOGGER.logMessage("Generating evolution for Q = %.9lf", Q);

		while((std::sqrt(t)*z) > std::sqrt(fac_cutoff * tMin)) {
			Emission emission = generateEmission(std::sqrt(t)*z, std::sqrt(tMin), fac_t, alphaSOver);
			t = emission.t;
			z = emission.z;
			Double pT_2 = emission.pT_2;
			Double m_2  = emission.m_2;
			Bool generated = emission.generated;
			Bool continueEvolution = emission.continueEvol;

			if (!continueEvolution) {
				LOGGER.logMessage("No further emissions.");
				LOGGER.logMessage("Total emissions: %u", numEmissions);
				LOGGER.logMessage("-------");
				return emissions;
			}

			if (t < (fac_cutoff*tMin)) {
				LOGGER.logMessage("Emission rejected at sqrt(t) = %.9lf since it is below cutoff", std::sqrt(t));
				LOGGER.logMessage("Total emissions: %u", numEmissions);
				return emissions;
			}

			if (z != 1.0) {
				emissions.emplace_back(std::sqrt(t), z, std::sqrt(pT_2), std::sqrt(m_2), true, true);
				numEmissions++;

				LOGGER.logMessage("Successful emission at sqrt(t) = %.9lf", std::sqrt(t));
				LOGGER.logMessage("  z = %.9lf", z);
				LOGGER.logMessage("  pT = %.9lf", std::sqrt(pT_2));
				LOGGER.logMessage("  m = %.9lf",  std::sqrt(m_2));
			}
		}

		LOGGER.logMessage("No futher emissions");
		LOGGER.logMessage("Total emission: %u", numEmissions);

		// simply list all of the emissions in a nice format
		logEmissionsTable(emissions);
		
		return emissions;
	}




	void PartonShower::logEmissionsTable(const std::vector<Emission>& emissions) {
		LOGGER.logMessage("Emissions table:");
		LOGGER.logMessage("+---+--------------------+---------------------+--------------------+---------------------------+");
		LOGGER.logMessage("| # |  Evo scale [GeV]   |         1-z         |      pT [GeV]      | virt. mass in a->bc [GeV] |");
		LOGGER.logMessage("+---+--------------------+---------------------+--------------------+---------------------------+");

		UInt32 count = 0;
		for (const Emission& e : emissions) {
			LOGGER.logMessage("| %d | %3d | %19.17lf | %18.16lf |     %17.15     |",
							  count++, e.t, 1.0-e.z, std::sqrt(e.pT_2), std::sqrt(e.m_2));
		}
		
		LOGGER.logMessage("+---+--------------------+---------------------+--------------------+---------------------------+");
	}
}; // namespace ColSim
