#ifndef __LHEF_HPP
#define __LHEF_HPP

#include "ColSim/Utils.hpp"

#include <utility>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>

namespace ColSim {
	namespace LHEF {

		// defines a ` name="xxxx" ` attribute
		typedef std::pair<std::string, std::string> XMLAttribute;


		// defines a single XML tag and its contents
		class LHEFTag {
		protected:
			// name and any attributes for this XLM tag
			std::string name;
			std::vector<XMLAttribute> attribs;

			// contents of tag as a string
			std::string contents;

			// any subtags within this main tag
			std::vector<LHEFTag*> subtags;

			// the i/o streams that the tag is read to/from
			std::ofstream outfile;
			std::ifstream infile;


		public:
			LHEFTag() {} // main constructor does nothing

			Optional<XMLAttribute> attrib(std::string const& name);
		};
	};
}; // namespace ColSim


#endif // __LHEF_HPP
