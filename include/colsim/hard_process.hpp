#ifndef __HARD_PROCESS_HPP
#define __HARD_PROCESS_HPP

#include "colsim/common.hpp"
#include "colsim/phase_space.hpp"
#include "colsim/particle.hpp"
#include "colsim/settings.hpp"

#include <limits>
#include <memory>

namespace colsim
{

	/** Struct to contain information related to the result
	 *  of the cross section calculation.
	 *  In particular, the result and error (obviously),
	 *  as well as maximum values of the weight and phase space points
	 *  to be used for hit-or-miss for the event generation.
	 */
	
	struct HardProcessResult
	{
		double result, error;
		double max_weight;
		std::vector<double> max_points;

		HardProcessResult() = delete;
		HardProcessResult(uint dims)
			: max_points(dims, 0.0)
		{}
		~HardProcessResult() = default;

		double dims() const { return max_points.size(); }
	};


	/** Base class containing information related to calculation
	 *  and generation of events pertaining to specific hard processes.
	 */
	class HardProcess
	{
	protected:
		/** The underlying phase space object tasked with the random
		 *  but self-consistent generation of phase space points
		 *  characteristic for this process
		 */
		std::unique_ptr<PhaseSpace> _phase_space;

	public:
		HardProcess() : _phase_space(nullptr) {}
		virtual ~HardProcess() = default;

		PhaseSpace const& get_phase_space() const { return *_phase_space; }
		PhaseSpace& get_phase_space() { return *_phase_space; }


		struct Result
		{
			double weight;
			std::vector<double> additional_vals;

			Result(double _weight, std::initializer_list<double> vals = {})
				: weight{_weight}, additional_vals{vals}
			{}


			static inline Result invalid_result()
			{
				return Result(std::numeric_limits<double>::max(), {});
			}

			bool operator==(Result const& other) const
			{
				if (weight == other.weight) {
					if (additional_vals == other.additional_vals)
						return true;
				}
				return false;
			}
		};

		
		/** Calculates the differential cross section
		 *  given an array of phase space points
		 */
		virtual Result dsigma(std::vector<double> const& phaseSpacePoints) = 0;
		
		/** Simulates the hard scattering process and calculates the
		 *  cross section and max weight values.
		 */
		HardProcessResult calculate();

		/** Generates (4) momenta and places them inside @a momenta. 
		 */
		virtual void generate_particles(std::vector<Particle>& momenta) = 0;
	};

	// p + p -> l + l
	class PP2Zg2ll : public HardProcess
	{
    public:
		PP2Zg2ll();
		~PP2Zg2ll() = default;

		Result dsigma(const std::vector<double>& phaseSpacePoints) override;

		void generate_particles(std::vector<Particle>& momenta) override;


	private:
		// helpers for calculation of the cross section
		double kappa() const;
		double chi1(double s_hat) const;
		double chi2(double s_hat) const;
		double a0(uint quarkType, double s_hat) const;
		double a1(uint quarkType, double s_hat) const;

		/** Computes the partonic cross section.
		 */
		double dsigma_hat(
			double cosTheta,
			uint quarkType,
			double s_hat);

		/** Computes the full weight (i.e. value of integrand).
		 *  In particular, this function calculates the PDF values.
		 */
		double compute_weight(
			double cosTheta,
			double x1, double x2,
			double s_hat);

	};


	class PP2Jets : public HardProcess
	{
	private:
		constexpr static double _eta4 = 0.5;
    public:
		PP2Jets();
		~PP2Jets() = default;

		Result dsigma(std::vector<double> const& phaseSpacePoints) override;
		void generate_particles(std::vector<Particle>& momenta) override;

	private:
		// some phase space calculations
		inline double calc_eta_star(double eta3, double eta4) {
			return (1.0/2.0)*(eta3 - eta4);
		}
		inline double calc_eta_bar(double eta3, double eta4) {
			return (1.0/2.0)*(eta3 + eta4);
		}

		// compute the momentum fractions from original phase space variables
		inline double calc_x1(double eta3, double eta4, double E_t) {
			return (2*E_t / SETTINGS.ecm) * std::cosh(calc_eta_star(eta3, eta4)) * std::exp(calc_eta_bar(eta3, eta4));
		}
		inline double calc_x2(double eta3, double eta4, double E_t) {
			return 2*E_t / SETTINGS.ecm * std::cosh(calc_eta_star(eta3, eta4)) * std::exp(-calc_eta_bar(eta3, eta4));
		}
		
		// s_hat, t_hat, u_hat calculations
		inline double calc_sh(double eta_star, double E_t) {
			return 4*E_t*E_t * std::cosh(eta_star)*std::cosh(eta_star);
		}
		inline double calc_th(double eta_star, double E_t) {
			return -2*E_t*E_t * std::cosh(eta_star)*std::exp(-eta_star);
		}
		inline double calc_uh(double eta_star, double E_t) {
			return -calc_sh(eta_star, E_t) - calc_th(eta_star, E_t);
		}

		// individual calculations of the partonic cross sections
		double qqprime2qqprime(double eta3, double eta4, double et, double alphas);
		double qqprimebar2qqprimebar(double eta3, double eta4, double et, double alphas);
		double qq2qq(double eta3, double eta4, double et, double alphas);
		double qqbar2qprimeqprimebar(double eta3, double eta4, double et, double alphas);
		double gg2gg(double eta3, double eta4, double et, double alphas);
		double qqbar2gg(double eta3, double eta4, double et, double alphas);
		double gg2qqbar(double eta3, double eta4, double et, double alphas);
		double qg2qg(double eta3, double eta4, double et, double alphas);
	};
	
}; // namespace ColSim



#endif // __HARD_PROCESS_HPP
