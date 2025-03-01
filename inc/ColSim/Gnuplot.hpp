#ifndef __GNUPLOT_HXX
#define __GNUPLOT_HXX

#include "ColSim/Types.hpp"

#include <vector>
#include <string>

namespace ColSim {

class Gnuplot {
		
private:
	static const std::string DEFAULT_SCRIPT_PATH;
	
	std::vector<std::string> commands; /*!< List of string commands before compilation into script. */
	UInt32 numBins;
	UInt32 min, max;


public:
	Gnuplot();
	~Gnuplot() {}

	/** Sets the title of the plot.
	 *  Can also be done via the constructor.
	 *
	 *  @param title Title to set for the plot.
	 */ 
	void setTitle(const std::string& title);

	void plotData(const std::string& datafile, UInt32 col);

	inline void setBinNum(const UInt32 _numBins) {
		numBins = _numBins;
	}
	inline void setMinMax(const UInt32 _min, const UInt32 _max) {
		min = _min;
		max = _max;
	}

	/** Saves all commands to a file and passes it to Gnuplot.
	 */ 
	void plot();

}; // class Gnuplot


} // namespace ColSim






#endif // endif
