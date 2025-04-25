#include "ColSim/Gnuplot.hpp"

#include "ColSim/Logger.hpp"
#include "ColSim/Utils.hpp"

#include <fstream>
#include <sstream>
#include <cstdlib>


namespace ColSim {


	Gnuplot::~Gnuplot() {
		// ensure the file streams are closed
		dataFile.close();
		scriptFile.close();
	}

	void Gnuplot::openDataFile(const std::string& fileName,
							   const std::vector<std::string>& colNames)
	{
		dataFile.open(fileName);
		if (!dataFile) {
			LOGGER.logAbort("File name %s for gnuplot data file cannot be opened.",
							fileName.c_str());
		}

		// top-level comment indicating column names
		dataFile << '#';
		for (const std::string& name : colNames)
			dataFile << name << '\t';
		dataFile << '\n';

		// store the column names to be used for generating plots/filenames
		this->colNames = colNames;
		dataFileName = fileName;
	}

	void Gnuplot::addDataPoint(const std::vector<Double>& point) {
		if (point.size() < 1) {
			LOGGER.logWarning("(Gnuplot::addDataPoint) found empty array of points");
			return;
		}
		// only add the point if it falls within the min/max
		bool valid = true;
		for (UInt i=0; i<point.size(); i++) {
			if ((point.at(i) < mins.at(i)) || (point.at(i) > maxes.at(i)))
				valid = false;
		}
		if (!valid)
			return;
		
		for (const Double p : point)
			dataFile << p << '\t';
		dataFile << '\n';
	}

	void Gnuplot::addDataPoints(const std::vector<std::vector<Double>>& points) {
		for (const std::vector<Double>& point : points)
			addDataPoint(point);
		dataFile << '\n'; // extra newline just in case

		// assume this is all of the data.
		dataFile.close();
	}



	void Gnuplot::plot() {
		Int res;
		
		// create a scripts directory
	    res = system("mkdir -p scripts plots");
		

		for (UInt i=0; i<colNames.size(); i++) {
			
			// do a copy of the name
			std::string name = colNames.at(i);
			
			// if parentheses in the name, replace with underscores
			// closing ones can be removed so there are no trailing underscores
			// also remove spaces
			// so, something like "y (rapidity)" becomes "y_rapidity"
			Replace(name, "(", "_");
			Replace(name, ")", "");
			Replace(name, " ", "");

			LOGGER.logMessage("Saving plot '%s'", name.c_str());

			std::string fileName = "scripts/" + name + ".gplt";
			scriptFile.open(fileName);
			if (!scriptFile)
				LOGGER.logAbort("File '%s' could not be opened.", fileName.c_str());
			
			// indicate output to a pdf
			scriptFile << "set terminal pdfcairo enhanced color notransparent\n";
			scriptFile << "set output 'plots/" << name << ".pdf'\n";

			// set title and axis labels
			// these are optional so check if they have been set yet
			if (titles.size() >= 1)
				scriptFile << "set title '" << titles[i] << "'\n";
			if (xlabels.size() >= 1)
				scriptFile << "set xlabel '" << xlabels[i] << "'\n";
			if (ylabels.size() >= 1)
				scriptFile << "set ylabel '" << ylabels[i] << "'\n";

			// set up histogram information
			const Double binWidth = deltas.at(i)/static_cast<Double>(numBins);
			scriptFile << "binwidth=" << binWidth << '\n';
			scriptFile << "hist(x,w)=w*floor(x/w)+w/2.0\n";
			scriptFile << "set boxwidth binwidth*0.9\n";
			scriptFile << "set style fill solid 0.5\n";


			// do the plotting
			scriptFile << "set xrange [" << mins.at(i) << ":" << maxes.at(i) << "]\n";
			scriptFile << "plot '" << dataFileName << "' using (hist($" << (i+1) << ",binwidth)):(1.0) "
					   << "smooth frequency with boxes lc rgb\"blue\" notitle\n";

			scriptFile.close();

			res = system(std::string("gnuplot " + fileName).c_str());
		}
		(void)res;
	}
	

} // namespace ColSim
