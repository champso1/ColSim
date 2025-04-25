#include <ColSim/ColSim.hpp>
#include <ColSim/Gnuplot.hpp>
using namespace ColSim;

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
using namespace std;

int main() {
	// required for setup of random number generator
	ColSimMain _;

	// create gnuplot object
	Gnuplot plot;

	// specify plot information
	vector<double> min{0.0}, max{3.14}, delta{3.14};
	double numBins = 25;
	plot.setHistInfo(min, max, delta, numBins);

	// give col/file name, open the datafile
	vector<string> colNames{"sin2_x"};
	plot.openDataFile("data.dat", colNames);


	// grab values of sin(x)
	// in general we take a vector of vectors
	// for multiple variables,
	// that's the reason for adding a vector
	uint numPoints = 0;
	while (numPoints < 100000) {
		double x = randDouble() * 3.14;
		double r = randDouble();
		double val = sin(x)*sin(x);
		double ratio =  val / 1.0;
		if (r < ratio) {
		    plot.addDataPoint(vector<double>{x});
			numPoints++;
		}
	}
	
	cout << "opened data file and specified column names\n";
	
	
	// define the x and y labels
	vector<string> xlabel{"x"}, ylabel{"sin^2(x)"}, title{"Histogram format for sin^2(x)"};
	plot.setTitles(title);
	plot.setXLabels(xlabel);
	plot.setYLabels(ylabel);

	cout << "added labels\n";

	// make the plots
	plot.plot();

	cout << "plotted\n";
	
   	return 0;
}
