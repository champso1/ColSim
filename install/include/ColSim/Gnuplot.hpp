#ifndef __GNUPLOT_HXX
#define __GNUPLOT_HXX

#include "ColSim/Types.hpp"

#include <vector>
#include <string>
#include <fstream>

namespace ColSim {

	class Gnuplot {
		
	private:
		// file streams for data and for the script
		std::ofstream dataFile, scriptFile;
		std::string dataFileName;

		// storage for the column names
		// to be used for generating the plots
		std::vector<std::string> colNames;

		// for determining histogram binning and ranges
		std::vector<Double> mins, maxes, deltas;
		UInt32 numBins;

		std::vector<std::string> titles, xlabels, ylabels;

	public:
		Gnuplot() { }
		~Gnuplot();

		void openDataFile(const std::string& fileName, const std::vector<std::string>& colNames);
		
		void addDataPoint(const std::vector<Double>& point);
		void addDataPoints(const std::vector<std::vector<Double>>& points);

		inline void setHistInfo(const std::vector<Double>& mins,
							    const std::vector<Double>& maxes,
								const std::vector<Double>& deltas,
								const UInt32 numBins)
		{
			this->mins = mins;
		    this->maxes = maxes;
			this->deltas = deltas;
			this->numBins = numBins;
		}


		inline void setTitles(const std::vector<std::string>& _titles) { titles = _titles; }
		inline void setXLabels(const std::vector<std::string>& _xlabels) { xlabels = _xlabels; }
		inline void setYLabels(const std::vector<std::string>& _ylabels) { ylabels = _ylabels; }

		
		/** Saves all commands to a file and passes it to Gnuplot.
		 *  Also flushes/closes file streams.
		 */ 
		void plot();

	}; // class Gnuplot


} // namespace ColSim






#endif // endif
