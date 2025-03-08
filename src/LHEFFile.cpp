#include "ColSim/LHEF.hpp"
#include "ColSim/Utils.hpp"
#include "ColSim/Settings.hpp"

#include <exception>
#include <stdexcept>
#include <sstream>


namespace ColSim {
	namespace LHEF {

		void LHEFile::open(const std::string& filePath, OpenMode mode) {
			// close both files first in the case that
			// one is already open
			outfileStream.close();
			infileStream.close();
			switch (mode) {
				case IN: {
					infileStream.open(filePath, std::ios_base::in);
					if (!infileStream)
						throw std::runtime_error("Could not open LHEFile\n");
				} break;
				case OUT:
					outfileStream.open(filePath, std::ios_base::out);
					if (!outfileStream)
						throw std::runtime_error("Could not open LHEFile\n");
					break;
			}
		}
		

		std::vector<Particle> LHEFile::getParticles() const {
			std::vector<Particle> particles;

			for (const Event& event : events) {
				for (const EventParticle& p : event.particles) {
					particles.emplace_back(
					    FourVector(p.e, p.px, p.py, p.pz),
					    p.id,
					    "X");
				}
			}

			return particles;
		}

		void LHEFile::readData() {
			if (!infileStream)
				throw std::runtime_error("Could not open LHEFile.\n");

			std::string line;
			Bool inEvent = false;
			UInt8 test = 0;

			Event event; // dummy event
			
			while (std::getline(infileStream, line)) {
				// if we find the start of en event,
				// flip the inEvent bit
				if (line.find("<event>") != std::string::npos) {
					inEvent = true;
					continue; // can go ahead and continue
				}

				// while we are in an event, 
				if (inEvent) {
					std::vector<std::string> tokens;
					SplitString(line, " ", tokens);

					switch (tokens.size()) {
						case 6:
							event.weight = std::stod(tokens[2]);
							break;
						case 13:
						    event.particles
								 .emplace_back(readParticle(tokens));
							break;
						case 4: {
							std::string key = tokens[1];
							Erase(tokens[1], "id=");
							Erase(key, ">");
							Erase(key, "'");
							event.multiWeights[key] = std::stod(tokens[2]);
						};
					}
				}


				// if we find the closing tag,
				// append to the list of events (i.e. copy)
				// then clear this local event to make room
				// for new events
				// and toggle the inEvent bit
				if (line.find("</event>") != std::string::npos) {
					events.emplace_back(event);
					event.clear();
					inEvent = false;

					// only read two events for testing
					if (test++ >= 2)
						break;
				}
			}
		}


		void LHEFile::writeData() {
			writeHeader();

			for (const Event& e : events) {
				std::vector<std::string> particleStrings;
				for (const EventParticle& p : e.particles) {
					std::stringstream particleString;
				    particleString << "TODO: particle :^)";
					particleStrings.emplace_back(particleString.str());
				}
				
				outfileStream << "<event>\n";
				for (const std::string& ps : particleStrings)
					outfileStream << "\t" << ps << "\n";
				outfileStream << "</event>\n";
			}
			
			writeFooter();
		}


		EventParticle LHEFile::readParticle(const std::vector<std::string>& line) {
			return EventParticle(std::stoi(line[0]), std::stoi(line[1]),
								 std::stod(line[6]), std::stod(line[7]),
								 std::stod(line[8]), std::stod(line[9]),
								 std::stod(line[10]));
		}

		void LHEFile::writeHeader() {
			const Double protonEnergy = SETTINGS.ECM/2.0;
			const Double xs = SETTINGS.crossSection;
			const Double stddev = SETTINGS.crossSectionError;
			outfileStream << "<LesHouchesEvents version=\"1.0\">\n";
			outfileStream << "<!--\n";
			outfileStream << "Generated with ColSim\n";
			outfileStream << "-->\n";
			outfileStream << "<init>\n";
			outfileStream << "\t11 -11 ";
			outfileStream << protonEnergy << " " << protonEnergy;
			outfileStream << "0 0 7 7 1 1\n";
			// TODO:
			// add cross section and standard deviation
			outfileStream << "\t" << xs << " " << stddev << " 1.00000 9999\n";
			outfileStream << "</init>\n";
		}

		void LHEFile::writeFooter() {
			outfileStream << "</LesHouchesFile>";
		}
		
	}; // namespace LHEF
}; // namespace ColSim
