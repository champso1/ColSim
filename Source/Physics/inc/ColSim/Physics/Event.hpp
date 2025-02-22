#ifndef __EVENT_HPP
#define __EVENT_HPP

#include "ColSim/Core.h"
#include "ColSim/Physics/Particle.hpp"

#include <vector>

namespace ColSim {
	struct Event {
		std::vector<Particle> particles;
	};
}; // namespace ColSim



#endif // __EVENT_HPP
