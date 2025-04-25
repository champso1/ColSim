#ifndef __EVENT_HPP
#define __EVENT_HPP

#include "ColSim/Constants.hpp"
#include "ColSim/Particle.hpp"

#include <iostream>
#include <initializer_list>
#include <vector>

namespace ColSim {
	class Event {
	private:
		Double weight;
		std::vector<Particle> particles;
		
	public:
		Event(Double w, const std::vector<Particle>& p)
			: weight(w), particles(p) {}
		Event(Double w, std::vector<Particle>&& p)
			: weight(w), particles(p) {}
		Event(Double w, std::initializer_list<Particle> p)
			: weight(w), particles(p) {}

		// default copy/move constructors
		Event(const Event& event) = default;
		Event(Event&& event) = default;
		
		friend std::ostream& operator<<(std::ostream& o, const Event& event) {
			o << "Weight: " << event.weight << "\n";
			for (const Particle& p : event.particles)
				o << p;
			return o;
		}

		inline const std::vector<Particle>& getParticles() const { return particles; }
	};
}; // namespace ColSim



#endif // __EVENT_HPP
