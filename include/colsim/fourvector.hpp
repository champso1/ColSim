#ifndef __FOUR_VECTOR_HPP
#define __FOUR_VECTOR_HPP

#include <cmath>
#include <array>
#include <initializer_list>
#include <format>
#include <sstream>

namespace colsim
{
	class FourVector
	{
	public:
		using value_type = std::array<double, 4>;

	private:
		value_type _data{};
		
	public:
		FourVector() = default;
		FourVector(value_type const& data) : _data{data} {}
		FourVector(std::initializer_list<double> il);

		inline value_type const& data() const { return _data; }

		inline double norm() const
		{
			double res = 0.0;
			for (double x : _data)
				res += x*x;
			return res;
		}

		inline double e()  const { return _data[0]; }
		inline double px() const { return _data[1]; }
		inline double py() const { return _data[2]; }
		inline double pz() const { return _data[3]; }

		inline double pt_2() const { return _data[1]*_data[1] + _data[2]*_data[2]; }
		inline double pt() const { return std::sqrt(pt_2()); }

		inline FourVector& zboost(double beta) {
			double gamma = std::sqrt(1.0F / (1.0F - pow(beta, 2)));
			// store these to avoid modifying while computing
			double e = _data[0];
			double pz = _data[3];

			_data[0] = gamma * e - gamma * beta * pz;
			_data[1] = -gamma * beta * e + gamma * pz;

			return *this;
		}
	};
}; // namespace colsim

namespace std
{
	template <>
	struct formatter<colsim::FourVector, char>
	{
		template<class ParseContext>
		constexpr ParseContext::iterator parse(ParseContext& ctx)
		{
			return ctx.end();
		}
	
		template<class FmtContext>
		FmtContext::iterator format(colsim::FourVector const& s, FmtContext& ctx) const
		{
			std::ostringstream out;
			auto data = s.data();
			out << std::format("[{}, {}, {}, {}]", data[0], data[1], data[2], data[3]);
			return std::ranges::copy(std::move(out).str(), ctx.out()).out;
		}
	};
}


#endif // __FOUR_VECTOR_HPP
