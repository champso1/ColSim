#ifndef __HARD_PROCESS_HPP
#define __HARD_PROCESS_HPP

#include "ColSim/Types.hpp"
#include "ColSim/PhaseSpace.hpp"

#include <memory>
#include <functional>

namespace ColSim {

	struct XSResult {
		Double result, error;
	};

	class HardProcess {
	public:
		std::unique_ptr<PhaseSpace> phaseSpace;

	public:
		HardProcess() : phaseSpace(nullptr) {}
		virtual ~HardProcess() = default;
		
		inline const Double * getPhaseSpaceDeltas() { return phaseSpace->getDeltas(); }
		const Double * getPhaseSpaceMins() { return phaseSpace->getMins(); }
		const Double * getPhaseSpaceMaxes() { return phaseSpace->getMaxes(); }
		virtual Double dSigma(Double phaseSpacePoints[]) = 0;

		// returns a lambda expression that evaluates the differential cross section
		// given the phase space points
		// allows Monte Carlo integration
		inline std::function<Double(Double[])> genDSigmaLambda() {
			return [this](Double* phaseSpacePoints) -> double {
				return this->dSigma(phaseSpacePoints);
			};
		}
	};

	// p + p -> l + l
	class PP2Zg2ll : public HardProcess {
    public:
		PP2Zg2ll();
		~PP2Zg2ll() = default;
		
		Double Kappa() const;
		Double Chi1(const Double s_hat) const;
		Double Chi2(const Double s_hat) const;
		Double A0(const UInt8 quarkType, const Double s_hat) const;
		Double A1(const UInt8 quarkType, const Double s_hat) const;

		Double dSigmaHat(const Double cosTheta, const UInt8 quarkType, const Double s_hat);
		Double computeWeight(const Double cosTheta, const Double x1, const Double x2, const Double s_hat);

		Double dSigma(Double phaseSpacePoints[]) override;
	};



	class ee2Zg2mumu : public HardProcess {
	public:
		ee2Zg2mumu();
		~ee2Zg2mumu() = default;

		Double dSigmaHat();
		Double computeWeight();
		Double dSigma(Double* phaseSpacePoints) override;
	};
};



#endif // __HARD_PROCESS_HPP
