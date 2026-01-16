#ifndef __UTILS_HPP
#define __UTILS_HPP

#include <array>
#include <string_view>
#include <format>
#include <string>
#include <iostream>

namespace colsim
{
	// colors for printing to the terminal
	inline constexpr char const* ANSI_COLOR_RED =         "\x1b[31m";
	inline constexpr char const* ANSI_COLOR_GREEN =       "\x1b[32m";
	inline constexpr char const* ANSI_COLOR_YELLOW =      "\x1b[33m";
	inline constexpr char const* ANSI_COLOR_BLUE =        "\x1b[34m";
	inline constexpr char const* ANSI_COLOR_MAGENTA =     "\x1b[35m";
	inline constexpr char const* ANSI_COLOR_CYAN =        "\x1b[36m";
	inline constexpr char const* ANSI_COLOR_RESET =       "\x1b[0m";
	
	enum LogType : int
	{
		LOG_DEBUG = 0,
		LOG_INFO,
		LOG_WARNING,
		LOG_ERROR,
		LOG_ERROR_NOQUIT,
		LOG_NUM_LOG_TYPES
	};

	inline std::array<std::string_view, LOG_NUM_LOG_TYPES> log_string_reps{"DEBUG", "INFO", "WARNING", "ERROR", "ERROR"};
	inline std::array<std::string_view, LOG_NUM_LOG_TYPES> log_string_colors{ANSI_COLOR_GREEN, ANSI_COLOR_RESET, ANSI_COLOR_YELLOW, ANSI_COLOR_RED, ANSI_COLOR_RED};

	// TODO: better handle setting global flags
	inline bool debug_flag{false};
	inline bool& getDebugFlag() { return debug_flag; }

	template <typename... TArgs>
	void log(int log_type, std::string_view prefix, std::format_string<TArgs...> fmt_string, TArgs&& ...args)
	{
		if (log_type == LOG_DEBUG && !debug_flag)
			return;
			
		std::string log_text = std::vformat(fmt_string.get(), std::make_format_args(args...));
		std::string all_text = std::format("[{}] {}: {}\n", log_string_reps[log_type], prefix, log_text);
		std::cout << all_text;

		if (log_type == LOG_ERROR)
			exit(EXIT_FAILURE);
	}
}

#endif // __UTILS_HPP
