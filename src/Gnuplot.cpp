#include "ColSim/Gnuplot.hpp"

#include "ColSim/Logger.hpp"

#include <fstream>
#include <sstream>
#include <cstdlib>


namespace ColSim {

const std::string Gnuplot::DEFAULT_SCRIPT_PATH = "./res/script.dem";
	

Gnuplot::Gnuplot() {
	commands.emplace_back("set terminal pdfcairo enhanced color notransparent\n");
	commands.emplace_back("set output 'res/output.pdf'\n");
}


void Gnuplot::setTitle(const std::string& title) {
	std::ostringstream commandStream;
	commandStream << "set title '" << title << "'\n";
	commands.emplace_back(commandStream.str());
}

void Gnuplot::plot() {
	std::ofstream fileStream;
	fileStream.open(DEFAULT_SCRIPT_PATH);
	// error out otherwise
	// log error does this already
	if (!fileStream) 
	    LOGGER.logAbort("Unable to open/create default gnuplot script file: %s\n", DEFAULT_SCRIPT_PATH.c_str());

	for (const std::string& str : commands) {
		fileStream << str;
	}

	fileStream.close();
	const std::string systemCommand = "gnuplot " + DEFAULT_SCRIPT_PATH;
	system(systemCommand.c_str());
}

} // namespace ColSim
