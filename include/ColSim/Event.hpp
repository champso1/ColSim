#ifndef __EVENT_HPP
#define __EVENT_HPP

#include "ColSim/Constants.hpp"
#include "ColSim/Particle.hpp"

#include <iostream>
#include <initializer_list>
#include <vector>
#include <cfloat>

namespace ColSim {
	class Event {
	private:
		Double weight;
		std::vector<Particle> particles;
		Double Q, cos_theta, y;
		
	public:
		Event(Double w, const std::vector<Particle>& p, Double _Q=0, Double _cos_theta=0, Double _y=0)
			: weight(w), particles(p), Q(_Q), cos_theta(_cos_theta), y(_y) {}
		Event(Double w, std::vector<Particle>&& p, Double _Q=0, Double _cos_theta=0, Double _y=0)
			: weight(w), particles(p), Q(_Q), cos_theta(_cos_theta), y(_y) {}
		Event(Double w, std::initializer_list<Particle> p, Double _Q=0, Double _cos_theta=0, Double _y=0)
			: weight(w), particles(p), Q(_Q), cos_theta(_cos_theta), y(_y) {}

		// default copy/move constructors
		Event(const Event& event) = default;
		Event(Event&& event) = default;

		inline Double getQ() const { return Q; }
		inline Double getCosTheta() const { return cos_theta; }
		inline Double getY() const { return y; }

		inline const std::vector<Particle>& getParticles() const { return particles; }
		
		friend std::ostream& operator<<(std::ostream& o, const Event& event) {
			o << "Weight: " << event.weight << "\n";
			for (const Particle& p : event.particles)
				o << p;
			return o;
		}


		static inline Event invalidEvent() {
			return Event(DBL_MAX, {}, DBL_MAX, DBL_MAX, DBL_MAX);
		}
	};
}; // namespace ColSim



#endif // __EVENT_HPP
