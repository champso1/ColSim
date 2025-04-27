#ifndef __HARD_PROCESS_HPP
#define __HARD_PROCESS_HPP

#include "ColSim/Types.hpp"
#include "ColSim/PhaseSpace.hpp"
#include "ColSim/Event.hpp"
#include "ColSim/Math.hpp"
#include "ColSim/Settings.hpp"

#include <memory>
#include <functional>

namespace ColSim {

	/** Struct to contain information related to the result
	 *  of the cross section calculation.
	 *  In particular, the result and error (obviously),
	 *  as well as maximum values of the weight and phase space points
	 *  to be used for hit-or-miss for the event generation.
	 */
	
	struct HardProcessResult {
		Double result, error;
		Double maxWeight;
		std::vector<Double> maxPoints;

		HardProcessResult() = delete;
		HardProcessResult(USize numDims)
			: maxPoints(numDims, 0.0)
		{}
		~HardProcessResult() = default;

		Double getDims() const { return maxPoints.size(); }
	};


	/** Base class containing information related to calculation
	 *  and generation of events pertaining to specific hard processes.
	 */
	class HardProcess {
	protected:
		/** The underlying phase space object tasked with the random
		 *  but self-consistent generation of phase space points
		 *  characteristic for this process
		 */
		std::unique_ptr<PhaseSpace> phaseSpace;

	public:
		HardProcess() : phaseSpace(nullptr) {}
		virtual ~HardProcess() = default;

		/** Getter for the underlying phase space.
		 *  Returns a const reference.
		 */
		const PhaseSpace& getPhaseSpace() const {
			return *phaseSpace;
		}

		/** Getter for the underlying phase space.
		 *  Returns a modifiable (non-const) reference.
		 */
		PhaseSpace& getPhaseSpace() {
			return *phaseSpace;
		}


		struct Result {
			Double weight;
			std::vector<Double> additionalVals;

			Result(const Double _weight, std::initializer_list<Double> vals = {})
				: weight(_weight), additionalVals(vals)
			{}


			static inline Result invalidResult() {
				return Result(DBL_MAX, {});
			}

			Bool operator==(const Result& other) {
				if (weight == other.weight) {
					if (additionalVals == other.additionalVals)
						return true;
				}
				return false;
			}
		};

		
		/** Calculates the differential cross section
		 *  given an array of phase space points
		 */
		virtual Result dSigma(const std::vector<Double>& phaseSpacePoints) = 0;

		
		/** Simulates the hard scattering process and calculates the
		 *  cross section and max weight values.
		 */
		HardProcessResult calculate();
		
		

		/** Generates (4) momenta and places them inside @a momenta. 
		 */
		virtual void generateParticles(std::vector<Particle>& momenta) = 0;
	};

	// p + p -> l + l
	class PP2Zg2ll : public HardProcess {
    public:
		PP2Zg2ll();
		~PP2Zg2ll() = default;

		Result dSigma(const std::vector<Double>& phaseSpacePoints) override;

		void generateParticles(std::vector<Particle>& momenta) override;


	private:
		// helpers for calculation of the cross section
		Double Kappa() const;
		Double Chi1(const Double s_hat) const;
		Double Chi2(const Double s_hat) const;
		Double A0(const UInt8 quarkType, const Double s_hat) const;
		Double A1(const UInt8 quarkType, const Double s_hat) const;

		/** Computes the partonic cross section.
		 */
		Double dSigmaHat(const Double cosTheta,
						 const UInt8 quarkType,
						 const Double s_hat);

		/** Computes the full weight (i.e. value of integrand).
		 *  In particular, this function calculates the PDF values.
		 */
		Double computeWeight(const Double cosTheta,
							 const Double x1, const Double x2,
							 const Double s_hat);

	};


	class PP2Jets : public HardProcess {
	private:
		constexpr const static Double _eta4 = 0.5;
    public:
		PP2Jets();
		~PP2Jets() = default;

		Result dSigma(const std::vector<Double>& phaseSpacePoints) override;

		void generateParticles(std::vector<Particle>& momenta) override;

	private:
		// some phase space calculations
		inline Double calcEtaStar(Double eta3, Double eta4) {
			return (1.0/2.0)*(eta3 - eta4);
		}
		inline Double calcEtaBar(Double eta3, Double eta4) {
			return (1.0/2.0)*(eta3 + eta4);
		}

		// compute the momentum fractions from original phase space variables
		inline Double calcX1(Double eta3, Double eta4, Double E_t) {
			return (2*E_t / SETTINGS.ECM) * cosh(calcEtaStar(eta3, eta4)) * exp(calcEtaBar(eta3, eta4));
		}
		inline Double calcX2(Double eta3, Double eta4, Double E_t) {
			return 2*E_t / SETTINGS.ECM * cosh(calcEtaStar(eta3, eta4)) * exp(-calcEtaBar(eta3, eta4));
		}
		
		// s_hat, t_hat, u_hat calculations
		inline Double calcSHat(Double eta_star, Double E_t) {
			return 4*E_t*E_t * cosh(eta_star)*cosh(eta_star);
		}
		inline Double calcTHat(Double eta_star, Double E_t) {
			return -2*E_t*E_t * cosh(eta_star)*exp(-eta_star);
		}
		inline Double calcUHat(Double eta_star, Double E_t) {
			return -calcSHat(eta_star, E_t) - calcTHat(eta_star, E_t);
		}



		// individual calculations of the partonic cross sections
		Double qqprime2qqprime(Double eta3, Double eta4, Double Et, Double alphaS);
		Double qqprimebar2qqprimebar(Double eta3, Double eta4, Double Et, Double alphaS);
		Double qq2qq(Double eta3, Double eta4, Double Et, Double alphaS);
		Double qqbar2qprimeqprimebar(Double eta3, Double eta4, Double Et, Double alphaS);
		Double gg2gg(Double eta3, Double eta4, Double Et, Double alphaS);
		Double qqbar2gg(Double eta3, Double eta4, Double Et, Double alphaS);
		Double gg2qqbar(Double eta3, Double eta4, Double Et, Double alphaS);
		Double qg2qg(Double eta3, Double eta4, Double Et, Double alphaS);
	};
	
}; // namespace ColSim



#endif // __HARD_PROCESS_HPP
