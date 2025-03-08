#ifndef __HARD_PROCESS_HPP
#define __HARD_PROCESS_HPP

#include "ColSim/Types.hpp"
#include "ColSim/PhaseSpace.hpp"
#include "ColSim/Event.hpp"
#include "ColSim/Math.hpp"

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
			: result(0.0), error(0.0), maxWeight(0.0), // probably unecessary...
			  maxPoints(numDims, 0.0)
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

		/** Calculates the differential cross section
		 *  given an array of phase space points
		 */
		virtual Double dSigma(const std::vector<Double>& phaseSpacePoints) = 0;
		
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

		Double dSigma(const std::vector<Double>& phaseSpacePoints) override;

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
	
}; // namespace ColSim



#endif // __HARD_PROCESS_HPP
