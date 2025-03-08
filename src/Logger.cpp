#include "ColSim/Logger.hpp"

#include <ctime>
#include <cstdlib>

#include <exception>

namespace ColSim {

	void Logger::initFile(const std::string& fp) {
		filePath = fp;
		// open the file stream
		outFileStream.open(filePath);
		if (!outFileStream)
			throw std::runtime_error("Unable to open log file " + filePath);
		useFile = true;
	}

	Logger::~Logger() {
		outFileStream.close();
	}



	void Logger::getLogPrefix(char* prefixBuf, USize prefixBufSize) {
		time_t t = time(nullptr);
		strftime(prefixBuf, prefixBufSize, "%D %T", gmtime(&t));
	}



	void Logger::logError(const char* fmt, ...) {
		// includes the time and the 
		char prefixBuf[64] = {0};
		getLogPrefix(prefixBuf, 64);

		char messageBuf[256] = {0};
		va_list args;
		va_start(args, fmt);
		vsnprintf(messageBuf, 256, fmt, args);
	    va_end(args);

		char buf[512] = {0};
		snprintf(buf, 512, ANSI_COLOR_RED "%s [ERROR] %s\n" ANSI_COLOR_RESET, prefixBuf, messageBuf);

		fprintf(stderr, "%s", buf);
		if (useFile)
			outFileStream << buf;
	}

	void Logger::logWarning(const char* fmt, ...) {
		char prefixBuf[64] = {0};
		getLogPrefix(prefixBuf, 64);

		char messageBuf[256] = {0};
		va_list args;
		va_start(args, fmt);
		vsnprintf(messageBuf, 256, fmt, args);
	    va_end(args);

		char buf[512] = {0};
		snprintf(buf, 512, ANSI_COLOR_YELLOW "%s [WARNING] %s\n" ANSI_COLOR_RESET, prefixBuf, messageBuf);

		fprintf(stderr, "%s", buf);
		if (useFile)
			outFileStream << buf;
	}

	void Logger::logMessage(const char* fmt, ...) {
		char prefixBuf[64] = {0};
		getLogPrefix(prefixBuf, 64);

		char messageBuf[256] = {0};
		va_list args;
		va_start(args, fmt);
		vsnprintf(messageBuf, 256, fmt, args);
	    va_end(args);

		char buf[512] = {0};
		snprintf(buf, 512, "%s [INFO] %s\n", prefixBuf, messageBuf);

		fprintf(stderr, "%s", buf);
		if (useFile)
			outFileStream << buf;
	}


	void Logger::logAbort(const char* fmt, ...) {
		char prefixBuf[64] = {0};
		getLogPrefix(prefixBuf, 64);

		char messageBuf[256] = {0};
		va_list args;
		va_start(args, fmt);
		vsnprintf(messageBuf, 256, fmt, args);
	    va_end(args);

		char buf[512] = {0};
		snprintf(buf, 512, ANSI_COLOR_RED "%s [ERROR] %s\n" ANSI_COLOR_RESET, prefixBuf, messageBuf);

		fprintf(stderr, "%s", buf);
		if (useFile)
			outFileStream << buf;

		exit(EXIT_FAILURE);
	}
	
}; // namespace ColSim
