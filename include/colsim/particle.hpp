#ifndef __PARTICLE_HPP
#define __PARTICLE_HPP

#include <format>
#include <string>

#include "colsim/fourvector.hpp"

namespace colsim
{
	class Particle
	{
	private:
		FourVector _momentum;
		int _pid;
		std::string _name;

	public:
		Particle(int pid, std::string const& name)
			: _pid{pid}, _name{name}
		{}
		Particle(FourVector const& momentum,
				 int pid, std::string const& name)
			: _momentum{momentum}, _pid{pid}, _name{name}
		{}

		inline FourVector const& momentum() const { return _momentum; }
		inline int pid() const { return _pid; }
		inline std::string const& name() const { return _name; }

		inline double e() const { return  _momentum.e(); }
		inline double px() const { return _momentum.px(); }
		inline double py() const { return _momentum.py(); }
		inline double pz() const { return _momentum.pz(); }
		
		inline double pt() const {
			return std::sqrt(px()*px() + py()*py());
		}
		
	};
}; // namespace colsim

namespace std
{
	template <>
	struct formatter<colsim::Particle, char>
	{
		template<class ParseContext>
		constexpr ParseContext::iterator parse(ParseContext& ctx)
		{
			return ctx.end();
		}
	
		template<class FmtContext>
		FmtContext::iterator format(colsim::Particle const& p, FmtContext& ctx) const
		{
			std::ostringstream out;
			out << std::format("{} ({}) {}", p.name(), p.pid(), p.momentum());
			return std::ranges::copy(std::move(out).str(), ctx.out()).out;
		}
	};
}

#endif // __PARTICLE_HPP
