#ifndef __LOGGER_HPP
#define __LOGGER_HPP

#include "ColSim/Types.hpp"

#include <string>
#include <fstream>
#include <cstdarg>

namespace ColSim {

	// terminal colors for pretty printing
#define ANSI_COLOR_RED         "\x1b[31m"
#define ANSI_COLOR_GREEN       "\x1b[32m"
#define ANSI_COLOR_YELLOW      "\x1b[33m"
#define ANSI_COLOR_BLUE        "\x1b[34m"
#define ANSI_COLOR_MAGENTA     "\x1b[35m"
#define ANSI_COLOR_CYAN        "\x1b[36m"

#define ANSI_COLOR_RESET       "\x1b[0m"

	

	// --------------------
	// Logger
	// --------------------
	class Logger {
	public:
		Logger() : useFile(false) {}
		~Logger();

		/** Singleton getter.
		 */
		static Logger& getInstance() {
			static Logger logger;
			return logger;
		}

		/** Since the logger is a singleton, this is
		 *  essentially the constructor.
		 */
		void initFile(const std::string& fp);

		/** Logs an error to stderr and the provided output file
		 *  Colors it in red (if your terminal supports colors)
		 *  Indicates something has gone seriously wrong,
		 *  but the program is able to continue.
		 */
		void logError(const char* fmt ...);

		/** Logs a warning to stderr and the provided output file
		 *  Colors it in yellow (if your terminal supports colors)
		 *  Indicates a minor problem has occurred but largely
		 *  doesn't affect the program itself; execution continues.
		 */
		void logWarning(const char* fmt ...);

		/** Logs a message to stderr and the provided output file
		 *  No color used.
		 *  Used to output information or debug messages.
		 *  Does not affect program execuation.
		 */
		void logMessage(const char* fmt ...);

		/** Works the same as @a logError, but this will
		 *  also halt program execution.
		 */
		void logAbort(const char* fmt, ...);
	private:
		// info related to file
		Bool useFile;
		std::string filePath;
		std::ofstream outFileStream;

		// get log message prefix consisting of the time
		// NOTE: assumes this buffer is allocated already
		void getLogPrefix(char* prefixBuf, USize prefixBufSize);
	};

	
#define LOGGER Logger::getInstance()
	
	
}; // namespace ColSim



#endif // __LOGGER_HPP
