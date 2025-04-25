#include "ColSim/Types.hpp"
#include "ColSim/Utils.hpp"

#include <cstring>
#include <cstdlib>
#include <iostream>
#include <algorithm>


namespace ColSim {
	void SplitString(const std::string& _str, const std::string& _delim, std::vector<std::string>& res) {
		const char* str = _str.c_str();
		const char* delim = _delim.c_str();

		// copy the string data into a mutable buffer
		// so that we can tokenize it
		char buf[256];
		memset(buf, 0, 256);
		strncpy(buf, str, (sizeof buf)-1);

		char* ptr;
		ptr = strtok(buf, delim);
		while (ptr != nullptr) {
			res.emplace_back(ptr);
			ptr = strtok(nullptr, delim);
		}
	}

	
	template <typename T>
	Bool FindInVec(const std::vector<T>& vec, T item) {
		if (std::find(vec.begin(), vec.end(), item) != vec.end())
			return true;
		return false;
	}


	template <typename T>
	UInt32 CountInVec(const std::vector<T>& vec, T item) {
		UInt32 count = 0;
		for (const T& x : vec)
			if (x == item)
				count++;

		return count;
	}


	std::string& Erase(std::string& str, const std::string& in) {
		UInt64 inPos = str.find(in);
		if ((inPos = std::string::npos))
			return str;
		str.erase(inPos, in.length());
		return str;
	}


	std::string& Replace(std::string& str,
						 const std::string& orig,
						 const std::string& rep)
	{
		USize origPos;
		while ((origPos = str.find(orig)) != std::string::npos) {
			str.replace(origPos, orig.length(), rep);
		}

		return str;
	}

	template <typename T>
	std::vector<T> Join(const std::vector<T>& v1, const std::vector<T>& v2) {
		std::vector<T> res;
		for (auto x : v1)
			res.emplace_back(x);
		for (auto x : v2)
			res.emplace_back(x);

		return res;
	}


	// template instantiation for a few specific types
	template Bool FindInVec<std::string>(const std::vector<std::string>&, std::string);
	template UInt32 CountInVec(const std::vector<std::string>&, std::string);
	template std::vector<Double> Join(const std::vector<Double>& v1, const std::vector<Double>& v2);
}
