#include "colsim/parton_shower.hpp"

#include <fstream>
#include <limits>

#include "colsim/alphas.hpp"
#include "colsim/math.hpp"
#include "colsim/utils.hpp"

namespace colsim
{

	
	double PartonShower::tgamma(double z, double alphas_over)
	{
		return -2.0*alphas_over*CF*std::log1p(-z);
	}

	double PartonShower::tgamma_inv(double r, double alphas_over)
	{
		return 1.0 - std::exp(-0.5*r / (CF*alphas_over));
	}

	double PartonShower::zp_over(double t)
	{
		double Qcut = SETTINGS.evol_energy_cutoff;
		return 1.0-std::sqrt(Qcut*Qcut/t);
	}

	double PartonShower::zm_over(double t)
	{
		double Qcut = SETTINGS.evol_energy_cutoff;
		return std::sqrt(Qcut*Qcut/t);
	}

	double PartonShower::Pqq(double z)
	{
		return CF*(1.0 + z*z) / (1.0-z);
	}
	
	double PartonShower::Pqq_over(double z)
	{
		return 2.0*CF / (1.0-z);
	}


	double PartonShower::emission_scale_func(double log_T_Q2, void* _params)
	{
		Params* params = reinterpret_cast<Params*>(_params);
		double Q = params->Q;
		double r = params->r;
		double alphas_over = params->alphas_over;

		double t = Q*Q * std::exp(log_T_Q2);
		double rho = tgamma(zp_over(t), alphas_over) - tgamma(zm_over(t), alphas_over);

		return log_T_Q2 - (1.0/rho)*std::log(r);
	}

    double PartonShower::get_t_emission(double Q, double Qcut, double r, double tfac, double alphas_over)
	{
		auto emissionFunc = [&](double x, void* y) -> double {
			return emission_scale_func(x, y);
		};
		Params params{Q, Qcut, r, alphas_over};

	    double min = std::log(tfac * Qcut*Qcut / (Q*Q));
		double max = 0.0f;
		double tolerance = 1e-3F;

	    double root = Bisection(emissionFunc, min, max, 1000, tolerance, reinterpret_cast<void*>(&params));

		double sol = Q*Q * std::exp(root);

		// if we get something bad
		if (std::abs(emission_scale_func(root, reinterpret_cast<void*>(&params))) > tolerance)
			return std::numeric_limits<double>::max();
		
		return sol;
	}



	double PartonShower::get_z_emission(double t, double r, double alphas_over)
	{
	    return tgamma_inv( 
			tgamma(zm_over(t), alphas_over) + r * (tgamma(zp_over(t), alphas_over) - tgamma(zm_over(t), alphas_over)),
			alphas_over);
	}
	
	double PartonShower::get_pt_2(double t, double z)
	{
		return z*z * std::pow(1.0-z,2)*t;
	}
	
	double PartonShower::get_m_2(double t, double z)
	{
		return z*(1.0-z)*t;
	}

	PartonShower::Emission PartonShower::generate_emission(double Q, double Qcut, double tfac, double alphas_over)
	{
		bool generated = true;
		double r1 = rand_double();
		double r2 = rand_double();
		double r3 = rand_double();
		double r4 = rand_double();

		double t = get_t_emission(Q, Qcut, r1, tfac, alphas_over);
		if (t == std::numeric_limits<double>::max()) 
			return Emission{t, 1.0, 0.0, 0.0, generated, false};
		

		double z = get_z_emission(t, r2, alphas_over);
		double pT_2 = get_pt_2(t, z);
		double m_2 = get_m_2(t, z);

		//LOGGER.logMessage("Candidate emission scale: sqrt(t) = %.9lf", std::sqrt(t));
		//LOGGER.logMessage("Candidate momentum fraction: z = %.9lf", z);
		//LOGGER.logMessage("Candidate transverse momentum squared: pt^2 = %.9lf", pT_2);
		
		if (pT_2 < 0.0) {
			generated = false;
			//LOGGER.logMessage("Emission rejected due to negative pT^2");
		}

		double splittingFnOverRatio = Pqq(z) / Pqq_over(z);
		if (splittingFnOverRatio < r3) {
		    //LOGGER.logMessage("Emission rejected due to splitting function overestimate: p = %.9lf R = %.9lf", splittingFnOverRatio, r3);
			generated = false;
		} else {
			//LOGGER.logMessage("Emission NOT rejected due to splitting function overestimate: p = %.9lf R = %.9lf", splittingFnOverRatio, r3);
		}

		double alphas_over_ratio = _alphas.alphas_actual(t,z) / alphas_over;
		if (alphas_over_ratio < r4) {
			generated = false;
			//LOGGER.logMessage("Emission rejected due to alpha_s overestimate: a_s = %.9lf R = %.9lf", alphaSOverRatio, r4);
		} else {
			//LOGGER.logMessage("Emission NOT rejected due to alpha_s overestimate: a_s = %.9lf R = %.9lf", alphaSOverRatio, r4);
		}

		if (!generated) {
			z = 1.0;
			pT_2 = 0.0;
			m_2 = 0.0;
		}

		return Emission{t, z, pT_2, m_2, generated, true};
	}


	// TODO: cleanup this function
	std::vector<PartonShower::Emission> PartonShower::evolve(double Q, double Qmin, double alphas_over) {
		double tmin = Qmin*Qmin;
		uint numEmissions = 0;

		std::vector<Emission> emissions;

		double fac_t = 3.9999;
		double fac_cutoff = 4.0;

		double t = Q*Q;
		double z = 1.0;

		//LOGGER.logMessage("Generating evolution for Q = %.9lf", Q);

		while((std::sqrt(t)*z) > std::sqrt(fac_cutoff * tmin)) {
			Emission emission = generate_emission(std::sqrt(t)*z, std::sqrt(tmin), fac_t, alphas_over);
			t = emission.t;
			z = emission.z;
			double pT_2 = emission.pT_2;
			double m_2  = emission.m_2;
			bool continueEvolution = emission.continue_evol;

			if (!continueEvolution) {
				// LOGGER.logMessage("No further emissions.");
				// LOGGER.logMessage("Total emissions: %u", numEmissions);
#ifdef DEBUG
				if (numEmissions > 0)
					logEmissionsTable(emissions);
#endif				
				return emissions;
			}

			if (t < (fac_cutoff*tmin)) {
				// LOGGER.logMessage("Emission rejected at sqrt(t) = %.9lf since it is below cutoff", std::sqrt(t));
				// LOGGER.logMessage("Total emissions: %u", numEmissions);
#ifdef DEBUG
				if (numEmissions > 0)
					logEmissionsTable(emissions);
#endif
				return emissions;
			}

			if (z < 0.0) {
				// LOGGER.logMessage("Emission rejected at sqrt(t) = %.9lf since z<0: z=%.9lf", std::sqrt(t), z);
			    // LOGGER.logMessage("Total emissions: %u", numEmissions);
#ifdef DEBUG
				if (numEmissions > 0)
					logEmissionsTable(emissions);
#endif
				return emissions;
			}

			if (z != 1.0) {
				emissions.emplace_back(std::sqrt(t), z, std::sqrt(pT_2), std::sqrt(m_2), true, true);
				numEmissions++;

				// LOGGER.logMessage("Successful emission at sqrt(t) = %.9lf", std::sqrt(t));
				// LOGGER.logMessage("  z = %.9lf", z);
				// LOGGER.logMessage("  pT = %.9lf", std::sqrt(pT_2));
				// LOGGER.logMessage("  m = %.9lf",  std::sqrt(m_2));
			}
		}

		// LOGGER.logMessage("No futher emissions");
		// LOGGER.logMessage("Total emission: %u", numEmissions);

		// simply list all of the emissions in a nice format
#ifdef DEBUG
		if (numEmissions > 0)
			logEmissionsTable(emissions);
#else
		(void)numEmissions;
#endif

		// LOGGER.logMessage("Evolution completed! Performed %d emissions.", emissions.size());
		return emissions;
	}




	void PartonShower::log_emissions_table(std::vector<Emission> const& emissions) {
		log(LOG_INFO, "PartonShower::log_emissions_table()", "Emissions table:");
		log(LOG_INFO, "PartonShower::log_emissions_table()", "+---+--------------------+---------------------+--------------------+---------------------------+");
		log(LOG_INFO, "PartonShower::log_emissions_table()", "| # |  Evo scale [GeV]   |         1-z         |      pT [GeV]      | virt. mass in a->bc [GeV] |");
		log(LOG_INFO, "PartonShower::log_emissions_table()", "+---+--------------------+---------------------+--------------------+---------------------------+");

		uint count = 0;
		for (Emission const& e : emissions) {
			log(LOG_INFO,
				"PartonShower::log_emissions_table()", "| {} |      {:8.3f}      | {:19.15f} | {:18.15f} |     {:17.14f}     |",
				count++, e.t, 1.0-e.z, std::sqrt(e.pT_2), std::sqrt(e.m_2));
		}
		
		log(LOG_INFO, "PartonShower::log_emissions_table()", "+---+--------------------+---------------------+--------------------+---------------------------+");
	}
}; // namespace ColSim
