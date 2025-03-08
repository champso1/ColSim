#ifndef __LHEF_HPP
#define __LHEF_HPP

#include "ColSim/Types.hpp"
#include "ColSim/Particle.hpp"

#include <ios>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <vector>
#include <fstream>
#include <unordered_map>

namespace ColSim {
	namespace LHEF {

		// contains information about a particle within an event
		struct EventParticle {
			Int32 id, status;
			Double px, py, pz;
			Double e, m;

			// not used yet
		    // Int32 helicity, color, relation;

			EventParticle() = default;
			EventParticle(UInt32 _id, Double _status,
						  Double _px, Double _py, Double _pz,
						  Double _e, Double _m)
				: id(_id), status(_status),
				  px(_px), py(_py), pz(_pz),
				  e(_e), m(_m)
			{}

			// simple print function 
			friend std::ostream& operator<<(std::ostream& o, const EventParticle& p) {
				o << "[" << p.id << ", " << p.status;
				o << std::scientific;
				o << p.e << ", " << p.px << ", " << p.py << ", " << p.pz << "] ";
				o << "(" << p.m << ")\n";

				// go back to not using scientific
				std::cout.unsetf(std::ios::floatfield);
				return o;
			}
		};

		struct Event {
			std::vector<EventParticle> particles;
			Double weight;
			std::unordered_map<std::string, Double> multiWeights;

			Event() {}
			Event(const std::vector<EventParticle> _particles,
				  const Double _weight,
				  const std::unordered_map<std::string, Double> _multiWeights)
				: particles(_particles), weight(_weight), multiWeights(_multiWeights)
			{}

			// explicit copy constructor for extra assurance
			Event(const Event& event)
				: particles(event.particles),
				  weight(event.weight),
				  multiWeights(event.multiWeights)
			{}

			// default destructor since all members
			// are trivially destructable
			~Event() = default;

			/** Sets
			 */
			inline void clear() {
				particles.clear();
				weight = 0.0;
				multiWeights.clear();
			}
		};





		class LHEFile {
		private:
			std::ofstream outfileStream;
			std::ifstream infileStream;

		public:
			std::vector<Event> events;

		public:
			LHEFile() {}
			~LHEFile() {
				outfileStream.close();
				infileStream.close();
			};


			enum OpenMode {
				IN,
				OUT
			};
			

			/** Opens an LHE file
			 *
			 *  @note Should not be compressed (e.g. via gzip or similar)
			 */
			void open(const std::string& filePath, OpenMode mode);

			/** Closes the current file streams.
			 */
			inline void close() {
				outfileStream.close();
				infileStream.close();
			}

			/** Return a vector of Particles
			 *  from the event record.
			 *  Particles are easier to work with
			 *  and have some nicer methods, so
			 *  this is very nice.
			 */
			std::vector<Particle> getParticles() const;

			/** Reads all events from the currently opened
			 *  LHE file.
			 */
			void readData();

			/** Adds a single event to the event record.
			 */
			inline void addEvent(const Event& event) {
				events.emplace_back(event);
			}

			/** Writes all stored events
			 *  to the currently opened LHE file.
			 */
			void writeData();
			
		private:
			/** Helper function to construct a Particle
			 *  given an array of string tokens
			 *  from a line in the LHE file.
			 */
			EventParticle readParticle(const std::vector<std::string>& line);

			/** Helper function to write all of the initial
			 *  headers and other fields into an LHE file
			 *  before all events are written.
			 */
			void writeHeader();

			/** Helper function to write any final
			 *  blocks to the LHE file
			 */
			void writeFooter();
		};

	}; // namespace LHEF
}; // namespace ColSim


#endif // __LHEF_HPP
