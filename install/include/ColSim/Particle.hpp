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

		inline double e() const { return momentum.e; }
		inline double px() const { return momentum.px; }
		inline double py() const { return momentum.py; }
		inline double pz() const { return momentum.pz; }
		
		inline double pT() const {
			return std::sqrt(px()*px() + py()*py());
		}
		
	};



}; // namespace ColSim

#endif // __PARTICLE_HPP
