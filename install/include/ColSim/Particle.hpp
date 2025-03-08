#ifndef __PARTICLE_HPP
#define __PARTICLE_HPP

#include "ColSim/Types.hpp"
#include "ColSim/FourVector.hpp"

#include <iostream>
#include <string>
#include <vector>


namespace ColSim {


	class Particle {
	private:
		FourVector momentum;
		Int32 pid;
		std::string name;


	public:
		Particle(const Int32 _pid, const std::string& _name)
			: pid(_pid), name(_name)
		{}
		Particle(const FourVector& _momentum,
				 const Int32 _pid, const std::string& _name)
			: momentum(_momentum), pid(_pid), name(_name)
		{}
		
		friend std::ostream& operator<<(std::ostream& o, const Particle& p) {
			o << p.name << "(" << p.pid << ") ";
			o << p.momentum;
			o << "\n";

			return o;
		}
		
		const static std::vector<std::string> ALL_PARTICLE_NAMES;
	};



}; // namespace ColSim

#endif // __PARTICLE_HPP
