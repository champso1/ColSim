#ifndef __UTILS_HPP
#define __UTILS_HPP

#include "ColSim/Types.hpp"

#include <cstdio>
#include <cstdlib>

#include <vector>
#include <string>



// terminal colors 
#define ANSI_COLOR_RED         "\x1b[31m"
#define ANSI_COLOR_GREEN       "\x1b[32m"
#define ANSI_COLOR_YELLOW      "\x1b[33m"
#define ANSI_COLOR_BLUE        "\x1b[34m"
#define ANSI_COLOR_MAGENTA     "\x1b[35m"
#define ANSI_COLOR_CYAN        "\x1b[36m"

#define ANSI_COLOR_RESET       "\x1b[0m"





#define ERROR(...)													\
	do {															\
		fprintf(stderr, ANSI_COLOR_RED"[ERROR] (%s in %d):\n",		\
			    __FUNCTION__, __LINE__);							\
		fprintf(stderr, __VA_ARGS__);								\
		fprintf(stderr, ANSI_COLOR_RESET);							\
		exit(EXIT_FAILURE);											\
	} while(0)

#define WARNING(...)												\
	do {															\
		fprintf(stderr, ANSI_COLOR_YELLOW"[WARNING] (%s in %d):\n",	\
			    __FUNCTION__, __LINE__);							\
		fprintf(stderr, __VA_ARGS__);								\
		fprintf(stderr, ANSI_COLOR_RESET);							\
	} while(0)

#define DEBUG(...)														\
	do {																\
		fprintf(stderr, ANSI_COLOR_CYAN"[DEBUG] (line %d in %s()):\n",	\
			    __LINE__, __FUNCTION__);								\
		fprintf(stderr, __VA_ARGS__);									\
		fprintf(stderr, ANSI_COLOR_RESET);								\
	} while(0)





#define PTR_ASSERT(ptr, ...)								\
	if (!(ptr)) {											\
	LOG_ERROR(__VA_ARGS__)									\
		exit(EXIT_FAILURE);									\
	}


#define ASSERT(stmnt, ...) PTR_ASSERT(stmnt, __VA_ARGS__)




namespace ColSim {
	/** Splits a string by a delimiter and stores the result in @a res.
	 * 
	 *  @note Does not care if @a res has other elements; it simply appends them.
	 */
	extern void SplitString(const std::string& str, const std::string& delim, std::vector<std::string>& res);

	/** Returns whether the specified element was found in the vector.
	 *
	 *  @note Assumes the type is equality-comparable
	 */
	template <typename T>
	extern Bool FindInVec(const std::vector<T>& vec, T item);

	/** Returns the number of instances the specified element was found in the vector.
	 *
	 *  @note Assumes the type is equality-comparable
	 */
    template <typename T>
	extern UInt32 CountInVec(const std::vector<T>& vec, T item);





	template <typename T>
	class Optional {
	private:
		T value;
		bool has_value;

	public:
		Optional(T const& val) : value(val), has_value(true) {}
		Optional() : has_value(false) {}

		explicit operator bool() const noexcept {
			return has_value;
		}

		Optional& operator=(T const& val) {
			value = val;
			has_value = true;
			return *this;
		}

		Bool hasVal() const noexcept { return has_value; }
		T val() const noexcept { return value; }
	};
}





#endif // __UTILS_HPP
