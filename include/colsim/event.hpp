#ifndef __EVENT_HPP
#define __EVENT_HPP

#include "colsim/particle.hpp"

#include <vector>
#include <format>

namespace colsim
{
	class Event
	{
	private:
		double _weight;
		std::vector<Particle> _particles{};
		double _Q, _cos_theta, _y;
		
	public:
		Event(double w, std::vector<Particle> const& p, double Q=0, double cos_theta=0, double y=0)
			: _weight{w}, _particles{p}, _Q{Q}, _cos_theta{cos_theta}, _y{y}
		{}
		constexpr Event(double w, double Q=0, double cos_theta=0, double y=0)
			: _weight{w}, _Q{Q}, _cos_theta{cos_theta}, _y{y}
		{}

		inline double weight() const { return _weight; }
		inline double Q() const { return _Q; }
		inline double cos_theta() const { return _cos_theta; }
		inline double y() const { return _y; }
		inline std::vector<Particle> const& particles() const { return _particles; }
	};

	// I am so lazy lmao
	template <double N>
	static constexpr Event _INVALID_EVENT_BASE{N, N, N, N};
	static constexpr Event INVALID_EVENT = _INVALID_EVENT_BASE<std::numeric_limits<double>::max()>;
}; // namespace ColSim

namespace std
{
	template<>
	struct formatter<colsim::Event, char>
	{
		bool quoted = false;
	
		template<class ParseContext>
		constexpr ParseContext::iterator parse(ParseContext& ctx)
		{
			return ctx.end();
		}
	
		template<class FmtContext>
		FmtContext::iterator format(colsim::Event const& e, FmtContext& ctx) const
		{
			std::ostringstream out;
			out << std::format("\nWeight: {}\n", e.weight());
			for (auto& p : e.particles())
				out << std::format("\t{}", p);
			return std::ranges::copy(std::move(out).str(), ctx.out()).out;
		}
	};
}

#endif // __EVENT_HPP
