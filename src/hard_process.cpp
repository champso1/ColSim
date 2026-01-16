#include "colsim/hard_process.hpp"

#include <cmath>
#include <memory>

#include "LHAPDF/LHAPDF.h"

#include "colsim/utils.hpp"
#include "colsim/settings.hpp"
#include "colsim/math.hpp"
#include "colsim/common.hpp"
#include "colsim/phase_space.hpp"

namespace colsim
{
	HardProcessResult HardProcess::calculate()
	{
		uint num_dims = _phase_space->dims();
		
		double weight_sum = 0.0F;
		double weight_squared_sum = 0.0F;

		HardProcessResult res(num_dims);
		
	    // setup vector for phase space points
		std::vector<double> points;

		// here we get delta values which we need for Monte Carlo
	    const std::vector<double>& deltas = _phase_space->deltas();
		
		for (int i=0; i<SETTINGS.num_iterations; i++) {
			_phase_space->fill_phase_space(points);
			
			Result dsigma_res = dsigma(points);

			// if it is invalid, we must redo
			if (dsigma_res == Result::invalid_result()) {
				i--;
				continue;
			}
			
			double weight = dsigma_res.weight;
    
			// multiply by the deltas of the independent variables
			for (uint j=0; j<num_dims; j++)
			    weight *= deltas[j];

			// add to total weight
			weight_sum += weight;
			weight_squared_sum += pow(weight, 2);

			// set maxes
			if (weight > res.max_weight) {
				res.max_weight = weight;
				for (uint j=0; j<num_dims; j++)
					res.max_points[j] = points[j];
			}
		}

		// doing divisions, so we want a Double
		double num_evals = static_cast<double>(SETTINGS.num_iterations);
		
		res.result = weight_sum/num_evals;
		double variance = weight_squared_sum/num_evals
			              - std::pow(weight_sum/num_evals,2);
		res.error = std::sqrt(variance/num_evals);

		// scale to picobarns
		res.result *= MAGIC_FACTOR;
		res.error *= MAGIC_FACTOR;
		
		return res;
	}



	

	PP2Zg2ll::PP2Zg2ll() {
		_phase_space = std::unique_ptr<PhaseSpace>(new PhaseSpace_TauYCosth());
	}
	
	double PP2Zg2ll::kappa() const {
		return std::sqrt(2.0)*FERMI_CONSTANT*Z_MASS_2 / (4.0*PI*ALPHA);
	}

	double PP2Zg2ll::chi1(const double s_hat) const {
		return kappa()*s_hat*(s_hat-Z_MASS_2) / (std::pow(s_hat-Z_MASS,2) + Z_WIDTH_2*Z_MASS_2);
	}
	
	double PP2Zg2ll::chi2(const double s_hat) const {
		return std::pow(kappa(),2) * std::pow(s_hat,2) / (std::pow(s_hat-Z_MASS_2,2) + Z_WIDTH_2*Z_MASS_2);
	}
	
	double PP2Zg2ll::a0(uint quarkType, double s_hat) const {
		double CAe = -0.5, CVe = -0.5 + 2.0*WEINBERG_ANGLE;
		double CVf, CAf, Qf;
		if (quarkType == 0) { // up-type
			CVf = 0.5 - (4.0/3.0)*WEINBERG_ANGLE;
			CAf = 0.5;
			Qf = 2.0/3.0;
		} else { // down-type
			CVf = -0.5 + (2.0/3.0)*WEINBERG_ANGLE;
			CAf = -0.5;
			Qf = -1.0/3.0;
		}

	    return Qf*Qf - 2.0*Qf*CVe*CVf*chi1(s_hat) + (CAe*CAe + CVe*CVe)*(CAf*CAf + CVf*CVf) * chi2(s_hat);
	}
	
	double PP2Zg2ll::a1(uint quarkType, double s_hat) const {
		double A_mu = -0.5, V_mu = -0.5 + 2.0*WEINBERG_ANGLE;
		double V_quark, A_quark, Q_quark;
		if (quarkType == 0) { // up-type
			V_quark = 0.5 -(4.0/3.0)*WEINBERG_ANGLE;
			A_quark = 0.5;
			Q_quark = 2.0/3.0;
		} else {
			V_quark = -0.5 + (2.0/3.0)*WEINBERG_ANGLE;
			A_quark = -0.5;
			Q_quark = -1.0/3.0;
		}

	    return -4.0*Q_quark*A_mu*A_quark*chi1(s_hat) + 8.0*A_mu*V_mu*A_quark*V_quark*chi2(s_hat);
	}


	double PP2Zg2ll::dsigma_hat(double cosTheta, uint quarkType, double s_hat) {
		return (2.0*PI*std::pow(ALPHA,2) / (4.0*3.0*s_hat)) * (a0(quarkType,s_hat)*(1.0 + std::pow(cosTheta,2)) + a1(quarkType,s_hat)*cosTheta);
	}

	double PP2Zg2ll::compute_weight(double s_hat, double x1, double x2, double cosTheta) {
		double weight = 0.0;
		LHAPDF::PDF& pdf = *SETTINGS.pdf;
		// up-type quarks
		weight += dsigma_hat(cosTheta , 0, s_hat) * ((pdf.xfxQ2(2 , x1, s_hat) * pdf.xfxQ2(-2, x2, s_hat)) + (pdf.xfxQ2(4 , x1, s_hat) * pdf.xfxQ2(-4, x2, s_hat)));	
		weight += dsigma_hat(-cosTheta, 0, s_hat) * ((pdf.xfxQ2(-2, x1, s_hat) * pdf.xfxQ2(2 , x2, s_hat)) + (pdf.xfxQ2(-4, x1, s_hat) * pdf.xfxQ2(4 , x2, s_hat)));
		// down-type quarks
		weight += dsigma_hat(cosTheta , 1, s_hat) * ((pdf.xfxQ2(1 , x1, s_hat) * pdf.xfxQ2(-1, x2, s_hat)) + (pdf.xfxQ2(3 , x1, s_hat) * pdf.xfxQ2(-3, x2, s_hat)));
		weight += dsigma_hat(-cosTheta, 1, s_hat) * ((pdf.xfxQ2(-1, x1, s_hat) * pdf.xfxQ2(1 , x2, s_hat)) + (pdf.xfxQ2(-3, x1, s_hat) * pdf.xfxQ2(3 , x2, s_hat)));

		return weight;
	}


	HardProcess::Result PP2Zg2ll::dsigma( std::vector<double> const& phaseSpacePoints) {
		double S = SETTINGS.s;
		double E_TR_2 = SETTINGS.trans_energy_2;

		// independent variables
		double cos_theta = phaseSpacePoints[0];
		double rho       = phaseSpacePoints[1];
		double rand_y    = phaseSpacePoints[2];

		// other variables
	    double jacobian = (E_TR_2) / (std::cos(rho)*std::cos(rho) * S);

		double s_hat = E_TR_2*std::tan(rho) + E_TR_2;

		double ymax   = -0.5*std::log(s_hat/S);
		double deltay = 2.0*ymax;
		double y      = (2.0*rand_y - 1.0)*ymax;

		double x1 = std::sqrt(s_hat/S)*std::exp(y);
		double x2 = std::sqrt(s_hat/S)*std::exp(-y);

		if ((x1 > 1.0 || x1 < 0.0) || (x2 > 1.0 || x2 <0.0)) {
			log(LOG_WARNING, "PP2Zg2ll::dsigma()", "Encountered invalid x1 or x2: ({:.4f},{:.4f}). Retrying...", x1, x2);
			return Result::invalid_result();
		}

		double weight = compute_weight(s_hat, x1, x2, cos_theta);
		// go ahead and scale it this way here
		weight *= (jacobian * deltay);
		weight /= (x1 * x2);

		// pass through the COM energy and the momentum fraction (x1)
		return HardProcess::Result(weight, {std::sqrt(s_hat), x1});
	}


	void PP2Zg2ll::generate_particles(std::vector<Particle>& particles) {
		double S = SETTINGS.s;
		double ECM = SETTINGS.ecm;
		double E_TR_2 = SETTINGS.trans_energy_2;

		std::vector<double> phaseSpacePoints;
		_phase_space->fill_phase_space(phaseSpacePoints);

		double cos_theta = phaseSpacePoints[0];
		double rho      = phaseSpacePoints[1];
		double rand_y   = phaseSpacePoints[2];

		// other variables
		double s_hat = E_TR_2*std::tan(rho) + E_TR_2;
		double Q = std::sqrt(s_hat);

		double ymax   = -0.5*std::log(s_hat/S);
		double y      = (2.0*rand_y - 1.0)*ymax;

		double x1 = std::sqrt(s_hat/S)*std::exp(y);
		double x2 = std::sqrt(s_hat/S)*std::exp(-y);

		// boost parameter
		double beta = (x2-x1)/(x2+x1);
		// more phase space points
		double phi = rand_double()*2.0*M_PI;
		double sinPhi = std::sin(phi);
		double cosPhi = std::cos(phi);
		double sinTheta = std::sqrt(1.0 - cos_theta*cos_theta);

		particles.emplace_back(
			FourVector{0.5*x1*ECM, 0.0, 0.0, 0.5*x1*ECM},
			1, "q1");
		particles.emplace_back(
			FourVector{0.5*x2*ECM, 0.0, 0.0, -0.5*x2*ECM},
			1, "q2");
		particles.emplace_back(
			FourVector{0.5*Q, 0.5*Q*sinTheta*cosPhi, 0.5*Q*sinTheta*sinPhi, 0.5*Q*cos_theta}.zboost(beta),
			13, "l1");
		particles.emplace_back(
			FourVector{0.5*Q, -0.5*Q*sinTheta*cosPhi, -0.5*Q*sinTheta*sinPhi, -0.5*Q*cos_theta}.zboost(beta),
			-13, "l2");
	}





	PP2Jets::PP2Jets() {
		_phase_space = std::unique_ptr<PhaseSpace>(new PhaseSpace_EtEta());
	}


	double PP2Jets::qqprime2qqprime(double eta3, double eta4, double E_t, double alpha_s) {
		double eta_star = calc_eta_star(eta3, eta4);
		double s = calc_sh(eta_star, E_t);
		double t = calc_th(eta_star, E_t);
		double u = calc_uh(eta_star, E_t);

		return (4*M_PI*alpha_s*alpha_s)/(9.0*s*s) * ((s*s +  u*u)/(t*t));
	}
	double PP2Jets::qqprimebar2qqprimebar(double eta3, double eta4, double E_t, double alpha_s) {
		return qqprime2qqprime(eta3, eta4, E_t, alpha_s);
	}

	// scattering of identical quarks
	double PP2Jets::qq2qq(double eta3, double eta4, double E_t, double alpha_s) {
		double eta_star = calc_eta_star(eta3, eta4);
		double s = calc_sh(eta_star, E_t);
		double t = calc_th(eta_star, E_t);
		double u = calc_uh(eta_star, E_t);

		return (4.0*M_PI*alpha_s*alpha_s)/(9.0*s*s) * ((s*s + u*u)/(t*t) + (t*t + u*u)/(s*s) - (2.0/3.0)*((s*s)/(u*t)));
	}

	double PP2Jets::qqbar2qprimeqprimebar(double eta3, double eta4, double E_t, double alpha_s) {
		double eta_star = calc_eta_star(eta3, eta4);
		double s = calc_sh(eta_star, E_t);
		double t = calc_th(eta_star, E_t);
		double u = calc_uh(eta_star, E_t);


		return (4.0*M_PI*alpha_s*alpha_s)/(9.0*s*s) * ((t*t + u*u)/(s*s));
	}


	double PP2Jets::gg2gg(double eta3, double eta4, double E_t, double alpha_s) {
		double eta_star = calc_eta_star(eta3, eta4);
		double s = calc_sh(eta_star, E_t);
		double t = calc_th(eta_star, E_t);
		double u = calc_uh(eta_star, E_t);

		return (9.0*M_PI*alpha_s*alpha_s) / (2.0 * s*s) * (3 - (t*u)/(s*s) - (s*u)/(t*t) - (s*t)/(u*u));
    
	}
	double PP2Jets::qqbar2gg(double eta3, double eta4, double E_t, double alpha_s) {
		double eta_star = calc_eta_star(eta3, eta4);
		double s = calc_sh(eta_star, E_t);
		double t = calc_th(eta_star, E_t);
		double u = calc_uh(eta_star, E_t);

		return (32.0*M_PI*alpha_s*alpha_s)/(27.0*s*s) * (u/t + t/u - (9.0/4.0)*((t*t + u*u)/(s*s)));
	}
	double PP2Jets::gg2qqbar(double eta3, double eta4, double E_t, double alpha_s) {
		double eta_star = calc_eta_star(eta3, eta4);
		double s = calc_sh(eta_star, E_t);
		double t = calc_th(eta_star, E_t);
		double u = calc_uh(eta_star, E_t);

		return (M_PI * alpha_s*alpha_s)/(6*s*s) * (u/t + t/u - (9.0/4.0)*((t*t + u*u)/(s*s)));
	}

	double PP2Jets::qg2qg(double eta3, double eta4, double E_t, double alpha_s) {
		double eta_star = calc_eta_star(eta3, eta4);
		double s = calc_sh(eta_star, E_t);
		double t = calc_th(eta_star, E_t);
		double u = calc_uh(eta_star, E_t);

		return (4.0*M_PI*alpha_s*alpha_s)/(9*s*s) * (-u/s - s/u + (9.0/4.0)*((s*s + u*u)/(t*t)));
	}


	HardProcess::Result PP2Jets::dsigma(std::vector<double> const& phaseSpacePoints) {
		double Et = phaseSpacePoints[0];
		double eta3 = phaseSpacePoints[1];

		// retrieve/calculate quantities required to grab PDF data
		double S = SETTINGS.s;
		LHAPDF::PDF& pdf = *SETTINGS.pdf;
		double alphaS = pdf.alphasQ2(S);
		double x1 = calc_x1(eta3, _eta4, Et);
		double x2 = calc_x2(eta3, _eta4, Et);
		const std::vector<uint> pids{1, 2, 3, 4, 5, 21}; // vector of PDG IDs for all quarks (but top) and gluon

		// must add all processes
		double total = 0.0;

		// q + q' -> q + q'
		double total_qqprime2qqprime = 0.0;
		for (uint i=0; i<pids.size(); i++) {
			for (int j=i+1; i<pids.size(); i++) {
				total_qqprime2qqprime += pdf.xfxQ2(pids[i], x1, S)*pdf.xfxQ2(pids[j], x2, S)
					* (1.0/M_PI)*qqprime2qqprime(eta3, _eta4, Et, alphaS);
			}
		}
		total += total_qqprime2qqprime;

		// q + q -> q + q
		double total_qq2qq = 0.0;
		for (uint i=0; i<pids.size(); i++) {
			total_qqprime2qqprime += pdf.xfxQ2(pids[i], x1, S)*pdf.xfxQ2(pids[i], x2, S)
				* (1.0/M_PI)*qq2qq(eta3, _eta4, Et, alphaS);
		}
		total += total_qq2qq;

		// q + qb -> q' + qb'
		double total_qqbar2qprimeqprimebar = 0.0;
		for (uint i=0; i<pids.size(); i++) {
			total_qqprime2qqprime += pdf.xfxQ2(pids[i], x1, S)*pdf.xfxQ2(pids[i], x2, S)
				* (1.0/M_PI)*qqbar2qprimeqprimebar(eta3, _eta4, Et, alphaS);
		}
		total += total_qqbar2qprimeqprimebar;

		
		// g + g -> g + g
	    total += pdf.xfxQ2(21, x1, S) * pdf.xfxQ2(21, x2, S) * gg2gg(eta3, _eta4, Et, alphaS);

		// q + qb -> g + g
		double total_qqbar2gg = 0.0;
		for (uint i=0; i<pids.size(); i++) {
			total_qqprime2qqprime += pdf.xfxQ2(pids[i], x1, S)*pdf.xfxQ2(pids[i], x2, S)
				* (1.0/M_PI)*qqbar2gg(eta3, _eta4, Et, alphaS);
		}
		total += total_qqbar2gg;

		// g + g -> q + qb
		double total_gg2qqbar = 0.0;
		for (uint i=0; i<pids.size(); i++) {
			total_qqprime2qqprime += pdf.xfxQ2(pids[i], x1, S)*pdf.xfxQ2(pids[i], x2, S)
				* (1.0/M_PI)*gg2qqbar(eta3, _eta4, Et, alphaS);
		}
		total += total_gg2qqbar;

		// q + g -> q + g
		double total_qg2qg = 0.0;
		for (uint i=0; i<pids.size(); i++) {
			total_qqprime2qqprime += pdf.xfxQ2(pids[i], x1, S)*pdf.xfxQ2(pids[i], x2, S)
				* (1.0/M_PI)*qg2qg(eta3, _eta4, Et, alphaS);
		}
		total += 2*total_qg2qg;

		total *= 2.0*Et;
		return Result(total);
	}

	void PP2Jets::generate_particles(std::vector<Particle>& momenta) {
		// does nothing yet
		(void)momenta;
		return;
	}

}; // namespace ColSim
