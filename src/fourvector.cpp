#include "colsim/fourvector.hpp"

#include <cstdlib>

namespace colsim
{
    FourVector::FourVector(std::initializer_list<double> il)
    {
        if (il.size() != 4)
            exit(EXIT_FAILURE);

        double const* data = std::data(il);
        for (uint i=0; i<_data.size(); ++i)
            _data[i] = data[i];
    }
}