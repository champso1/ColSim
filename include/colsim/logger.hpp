#pragma once

#include "colsim/common.hpp"

#include <format>
#include <string_view>
#include <array>
#include <iostream>
#include <unordered_map>
#include <mutex>

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
	inline constexpr char const* ANSI_LINEFEED_CLEAR =    "\033[2K";
	inline constexpr auto ANSI_LINEFEED_UP   = [](uint count){ return std::format("\033[{}F", count); };
	inline constexpr auto ANSI_LINEFEED_DOWN = [](uint count){ return std::format("\033[{}E", count); };
	inline constexpr std::string_view loading_block("█");

	
	/** Enum for defining a set of standard logging types. */
	enum LogType : int
	{
		LOG_NONE = -1,
		LOG_DEBUG,
		LOG_INFO,
		LOG_WARNING,
		LOG_ERROR,
		LOG_ERROR_NOQUIT,
		LOG_THREAD,
		LOG_NUM_LOG_TYPES
	};
	inline std::array<std::string_view, LOG_NUM_LOG_TYPES> log_string_reps{
		"DEBUG",
		"INFO",
		"WARNING",
		"ERROR",
		"ERROR",
		"THREAD"};
	inline std::array<std::string_view, LOG_NUM_LOG_TYPES> log_string_colors{
		ANSI_COLOR_GREEN,
		ANSI_COLOR_RESET,
		ANSI_COLOR_YELLOW,
		ANSI_COLOR_RED,
		ANSI_COLOR_RED,
		ANSI_COLOR_CYAN};

	/** @brief struct to store flags/options for logging */
	struct LogOptions final
	{
		int verbosity{LOG_INFO}; //!< anything equal to or above (in terms of priority) is displayed, everything else isn't
		bool use_log_output_stream{false}; //!< switch for whether we are logging to another output stream, like a file
		std::reference_wrapper<std::ostream> log_output_stream{std::ref(std::cerr)}; //!< actual output stream. only used if @a use_log_output_stream is true
	};
	/** @brief global options struct */
	inline LogOptions _log_options{};
	/** @brief returns the global options for setting values */
	inline LogOptions& getLogOptions() { return _log_options; }

	/** @brief Exception thrown for unrecoverable errors, i.e. every message logged with @a LOG_ERROR */
	class Error : public std::runtime_error
	{
	public:
		using std::runtime_error::runtime_error;
	};

	/**
	 *  @brief Prints a message to standard out and possibly an additional stream with a nice prefix.
	 *  @param log_type a type as in the enum @a LogType
	 *  @param prefix an identifier of some kind to include within the message, often a function name
	 *  @param fmt_string the format string of the log message, as used in std::format
	 *  @param args args used to format the string, as used in std::format
	 *  @throws Error if @a log_type is @a LOG_ERROR (regardless of the verbosity)
	 */
	template <typename... TArgs>
	void log(int log_type, std::string_view prefix, std::format_string<TArgs...> fmt_string, TArgs&& ...args)
	{
		int verbosity = getLogOptions().verbosity;
		bool shown = verbosity != LOG_NONE
			// we treat LOG_THREAD identical to LOG_INFO so we must handle that
			&& !(log_type == LOG_THREAD && verbosity > LOG_INFO)
			&& log_type >= verbosity;
		if (!shown && log_type != LOG_ERROR)
			return;

		std::string log_text = std::vformat(fmt_string.get(), std::make_format_args(args...));
		if (shown) {
			std::string all_text = std::format("{}[{}] {}: {}{}\n",
				log_string_colors[log_type], log_string_reps[log_type], prefix, log_text, ANSI_COLOR_RESET);
			if (getLogOptions().use_log_output_stream)
				getLogOptions().log_output_stream.get() << all_text;
			std::cout << all_text;
		}

		if (log_type == LOG_ERROR)
			throw Error(std::format("{}: {}", prefix, log_text));
	}

	/**
	 *  @brief A simple wrapper for printing to standard out and possibly an additional stream.
	 *  @param fmt_string the format string of the log message, as used in std::format
	 *  @param args args used to format the string, as used in std::format
	 */
	template <typename... TArgs>
	void log(std::format_string<TArgs...> fmt_string, TArgs&&... args)
	{
		if (getLogOptions().verbosity == LOG_NONE || getLogOptions().verbosity > LOG_INFO)
			return;
		
		std::string log_text = std::vformat(fmt_string.get(), std::make_format_args(args...));
		if (getLogOptions().use_log_output_stream)
			getLogOptions().log_output_stream.get() << log_text;
		std::cout << log_text;
	}

	inline std::unordered_map<uint, uint> log_threads_line_offset{};
	inline std::mutex log_threads_mutex; //!< guards @a log_threads_line_offset and the threaded progress output
	inline void registerThreadLogs(std::vector<uint> const& ids)
	{
		if (getLogOptions().verbosity == LOG_NONE || getLogOptions().verbosity > LOG_INFO)
			return;

		std::lock_guard<std::mutex> guard{log_threads_mutex};
		for (uint i=0; i<ids.size(); ++i) {
			log_threads_line_offset[ids[i]] = ids.size()-i;
			std::cout << '\n';
		}
	}
	/**
	 *  @brief Prints messages to standard out for threaded iterations in a nice way that doesn't flood stdout
	 *  @param thread_idx used in threaded printing, 
	 *  @param val the i in i/j, for the iteration count
	 *  @param end the j in i/j, for the iteration count
	 *  @param prefix an identifier of some kind to include within the message, often a function name
	 */
	template <typename... TArgs>
	void logThreadIterations(uint thread_idx, uint val, uint end, std::string_view prefix)
	{
		if (getLogOptions().verbosity == LOG_NONE || getLogOptions().verbosity > LOG_INFO)
			return;
		
		auto log_type = LOG_THREAD;

		uint count{};
		{
			std::lock_guard<std::mutex> guard{log_threads_mutex};
			if (auto it = log_threads_line_offset.find(thread_idx); it != log_threads_line_offset.end())
				count = it->second;
		}
		auto ansi_jump_up   = ANSI_LINEFEED_UP(count);
		auto ansi_jump_down = ANSI_LINEFEED_DOWN(count);
		double ratio = static_cast<double>(val)/static_cast<double>(end);
		int num_blocks = static_cast<int>(ratio*50.0);
		std::string blocks{};
		for (int i=0; i<num_blocks; ++i)
			blocks += std::string(loading_block);

		std::string all_text = std::format("{}[{}] {}: [{}] Iteration {:0>2}/{} ({: >3}%) [{: <50}]",
			log_string_colors[log_type],
			log_string_reps[log_type], prefix,
			thread_idx, val, end, static_cast<uint>(ratio*100.0),
			blocks,
			ANSI_COLOR_RESET);
		if (getLogOptions().use_log_output_stream)
			getLogOptions().log_output_stream.get() << all_text;

		{
			std::lock_guard<std::mutex> guard{log_threads_mutex};

			std::cout << ansi_jump_up << ANSI_LINEFEED_CLEAR
					  << all_text << std::flush
					  << ansi_jump_down << std::flush;
		}
	}

	inline void startLogIterations(){ return; }
	inline void endLogIterations()
	{
		if (getLogOptions().verbosity == LOG_NONE || getLogOptions().verbosity > LOG_INFO)
			return;
		
		std::cout << '\n';
	}
	/**
	 *  @brief Prints messages to standard out for iterations in a nice way that doesn't flood stdout
	 *  @param val the i in i/j, for the iteration count
	 *  @param end the j in i/j, for the iteration count
	 *  @param prefix an identifier of some kind to include within the message, often a function name
	 */
	template <typename... TArgs>
	void logIterations(uint val, uint end, std::string_view prefix)
	{
	    if (getLogOptions().verbosity == LOG_NONE || getLogOptions().verbosity > LOG_INFO)
			return;
		
		auto log_type = LOG_INFO;

		double ratio = static_cast<double>(val)/static_cast<double>(end);
		int num_blocks = static_cast<int>(ratio*50.0);
		std::string blocks{};
		for (int i=0; i<num_blocks; ++i)
			blocks += loading_block;

		std::string all_text = std::format("\r{}[{}] {}: Iteration {:0>2}/{} ({: >3}%) [{: <50}]{}",
			log_string_colors[log_type],
			log_string_reps[log_type], prefix,
			val, end, static_cast<uint>(ratio*100.0),
			blocks,
			ANSI_COLOR_RESET);
		if (getLogOptions().use_log_output_stream)
			getLogOptions().log_output_stream.get() << all_text;
		std::cout << all_text << std::flush;
	}
}
