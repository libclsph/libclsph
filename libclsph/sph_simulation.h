#ifndef _SPH_SIMULATION_H_
#define _SPH_SIMULATION_H_

#include <functional>
#include "common/structures.h"
#include "scene.h"
#include "../util/profile.hpp"

class sph_simulation {
public:
	sph_simulation() :
		parameters(),
		write_intermediate_frames(false) { }

	void simulate(int frame_count = 0);

	simulation_parameters parameters;
	precomputed_kernel_values precomputed_terms;

	std::function<void(particle*, const simulation_parameters&, bool, profile_data&)> pre_frame;
	std::function<void(particle*, const simulation_parameters&, bool, profile_data&)> post_frame;

	void load_settings(std::string fluid_file_name, std::string parameters_file_name);

	bool write_intermediate_frames;
    bool serialize;
	float initial_volume;
	scene current_scene;

private:
    void init_particles(particle* buffer,const simulation_parameters&);

	void sort_particles(
		particle*,
	    cl::Buffer&, cl::Buffer&,
	    cl::Kernel&, cl::Kernel&,
	    cl::CommandQueue&, cl::Context,
	    unsigned int*);

	void simulate_single_frame(
		particle*, particle*,
	    cl::Buffer&, cl::Buffer&,
        cl::Kernel&, cl::Kernel&, cl::Kernel&, cl::Kernel&,
	    cl::Kernel&, cl::CommandQueue&, cl::Context);

	profile_data profile;
};

#endif
