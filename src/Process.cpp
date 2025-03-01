namespace ColSim {

	/*
	// NOTE: this function is an absolute disaster at the moment
	// once more functions are implemented it will make more sense
	Process::Process(const std::string& process, const std::string& _config_file)
		: config_file(_config_file)
	{
		// determine the particles in the initial and final states
		std::vector<std::string> tokens;
		SplitString(process, " ", tokens);

		// ensure all particles are valid particles
		Bool found = true;
		for (const std::string& particle : tokens) {
			if (particle == "->") continue; // skip the arrow sign
			if (!FindInVec(Particle::ALL_PARTICLE_NAMES, particle)) {
					found = false;
					break;
				}
		}
		if (found)
			DEBUG("All particles are valid.\n");

		// ensure only 2->2 processes
		auto arrow_pos = std::find(tokens.begin(), tokens.end(), "->");
		USize num_initial_state_particles = std::distance(tokens.begin(), arrow_pos);
		USize num_final_state_particles = std::distance(arrow_pos, tokens.end());

		// only accept two particle initial/final states
		if ((num_initial_state_particles != 2) || (num_final_state_particles != 2))
			throw std::invalid_argument("We don't support anything other than 2->2 processes.");

	    // only accept pp collisions
		if ((tokens[0] != "p") || (tokens[1] != "p"))
		    throw std::invalid_argument("We don't support non-pp collisions (yet)");

		
		partonicXSFunc = PPtoDilepton;
		pdf = LHAPDF::mkPDF("cteq6l1", 0);
	}


	IntegrationResult Process::integrate() {
		// note that third one, y, has more scaling in the function
		// so we basically just make it 0-1 here
		// so that the monte carlo integration routine
		// doesn't do any scaling itself
		Double min[] = {
			-1.0F, RHO_MIN, 0.0F
		};
		Double max[] = {
			1.0F, RHO_MAX, 1.0F
		};
		MonteCarloParams params;
		params.num_evals = 1000000;
		params.num_dims = 3;
		params.min = min;
		params.max = max;
		params.extra_params = (void*)pdf;
		params.func = partonicXSFunc;
		
		IntegrationResult result = MonteCarloIntegrate(params);

		// scale to picobarns
		result.result *= MAGIC_FACTOR;
		result.variance *= MAGIC_FACTOR;
		result.error *= MAGIC_FACTOR;
		
		return result;
	}
	*/

	
}; // namespace ColSim
