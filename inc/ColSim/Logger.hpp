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

		static Logger& getInstance() {
			static Logger logger;
			return logger;
		}

		// set a file to be used
		void initFile(const std::string& filePath = DEFAULT_PATH);

		// message logging, optionally to the provided file
		void logError(const char* fmt ...);
		void logWarning(const char* fmt ...);
		void logMessage(const char* fmt ...);

		// also stops program execution
		void logAbort(const char* fmt, ...);
	private:
		// default path for logging to a file
		// check Logger.cpp for actual value
		static const std::string DEFAULT_PATH;

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
