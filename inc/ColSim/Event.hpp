#ifndef __EVENT_HPP
#define __EVENT_HPP

#include "ColSim/Particle.hpp"

#include <vector>

namespace ColSim {
	struct Event {
		std::vector<Particle> particles;
	};
}; // namespace ColSim



#endif // __EVENT_HPP
