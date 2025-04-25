#include <ColSim/ColSim.hpp>
using namespace ColSim;

#include <string>
using namespace std;

int main() {

	LOGGER.logMessage("This is a message.");
	LOGGER.logWarning("This is a warning.");
	LOGGER.logError("This is an error that will NOT abort the program.");

	string str = "this is a string!";
	LOGGER.logMessage("String: %s", str.c_str());

	int x = 42;
	LOGGER.logMessage("Integer: %d", x);

	double y = 5.679122;
	LOGGER.logMessage("Double: %.5lf", y);

	LOGGER.logAbort("This is an error that will abort the program.");
   
	
   	return 0;
}
