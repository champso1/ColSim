#define DEBUG

#include "ColSim/ColSim.hpp"
#include <ColSim/PartonShower.hpp>
using namespace ColSim;
using namespace std;

#include <cmath>
#include <string>
#include <fstream>
#include <vector>

int main() {
	ColSimMain colsim;
	colsim.init(ColSim::ColSimMain::PARTON_SHOWERING);

	vector<vector<double>> pts;
	
	vector<double> cutoffEnergies{1.0, 5.0, 10.0};
	
	for (const double e : cutoffEnergies) {
		SETTINGS.readString("EvolutionEnergyCutoff=" + to_string(e));
	
		colsim.start();
		colsim.generateEvents(10000);

		const vector<vector<PartonShower::Emission>>& emissionRecord = colsim.getEmissionRecord();
		vector<double> pt;

	    for (const vector<PartonShower::Emission>& emissions : emissionRecord) {
			for (const PartonShower::Emission& emission : emissions)
				pt.push_back(sqrt(emission.pT_2));
		}
		pts.emplace_back(pt);
	}
	colsim.stop();

	ofstream outfile("data-partonshower.dat");
	for (uint i=0; i<pts.at(0).size(); i++) {
		outfile << pts.at(0).at(i) << '\t';
			
		if (i < pts.at(1).size())
			outfile << pts.at(1).at(i) << '\t';
		if (i < pts.at(2).size())
			outfile << pts.at(2).at(i) << '\t';

		outfile << '\n';
	}
    
   	return 0;
}
