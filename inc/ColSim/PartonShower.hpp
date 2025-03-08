#ifndef __PARTON_SHOWER_HPP
#define __PARTON_SHOWER_HPP

#include "ColSim/Types.hpp"
#include "ColSim/Constants.hpp"

#include <random>

namespace ColSim {



	class PartonShower {
	protected:
		AlphaS alphaS;

	public:
		PartonShower() : alphaS(1) {};
		~PartonShower() = default;

		struct Emission {
			Double t, z, pT_2, m_2;
			Bool generated, continueEvol;

			Emission() = default;
			Emission(Double _t, Double _z, Double _pT_2, Double _m_2,
					 Bool _gen, Bool _cont)
				: t(_t), z(_z), pT_2(_pT_2), m_2(_m_2),
				  generated(_gen), continueEvol(_cont)
			{}
		};


		struct Params {
			Double Q, Qcut, r, alphaSOver;

			Params() = default;
			Params(Double _Q, Double _Qcut, Double _r, Double _alphaSOver)
				: Q(_Q), Qcut(_Qcut), r(_r), alphaSOver(_alphaSOver)
			{}
		};

	protected:
		Double tGamma(Double z, Double alphaSOver);
		Double tGammaInverse(Double r, Double alphaSOver);

		Double zpOver(Double t, Double Qcut);
		Double zmOver(Double t, Double Qcut);

		Double Pqq(Double z);
		Double PqqOver(Double z);

		Double emissionScaleFunc(Double log_T_Q2, void* _params);
		Double getTEmission(Double Q, Double Qcut, Double r, Double tfac, Double alphaSOver);


		Double getZEmission(Double t, Double Qcut, Double r, Double alphaSOver);
		Double getPT_2(Double t, Double z);
		Double getM_2(Double t, Double z);


		Emission generateEmission(Double Q, Double Qcut, Double tfac, Double alphaSOver);


	public:
		std::vector<Emission> Evolve(Double Q, Double Qmin, Double alphaSOver);
		


	private:
		// helper for printing the list of emissions
		// for a particle evolution
		// uses the logger's output stream
		void logEmissionsTable(const std::vector<Emission>& emissions);
		
	};

	class GluonShower : public PartonShower { };
	class PhotonShower : public PartonShower { };
	
}; // namespace ColSim

#endif // __PARTON_SHOWER_HPP
