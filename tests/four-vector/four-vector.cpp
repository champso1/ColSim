#include <ColSim/FourVector.hpp>
using namespace ColSim;

#include <iostream>
using namespace std;

int main() {

	FourVector v(89.0, 5.9, 3.4, 9.8);
	cout << v << '\n';
	cout << v.norm() << '\n';;
	cout << v.pT2() << '\n';
	cout << v.zBoost(0.5) << '\n';
	
   	return 0;
}
