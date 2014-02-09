#ifndef _SPH_SIMULATION_H_
#define _SPH_SIMULATION_H_

#include <functional>

#include "common/structures.h"

class sph_simulation {
public:
	sph_simulation() :
		parameters(),
		volumes(),
		write_intermediate_frames(false) { }

	void simulate(int frame_count = 0);

	simulation_parameters parameters;
	collision_volumes volumes;

	std::function<void(particle*, const simulation_parameters&, bool)> pre_frame;
	std::function<void(particle*, const simulation_parameters&, bool)> post_frame;

	void load_scene(std::string scene_file_name);
	void load_settings(std::string fluid_file_name, std::string parameters_file_name);

	bool write_intermediate_frames;

private:
	void init_particles(particle* buffer, int count);

	void sort_particles(
		particle*,
	    cl::Buffer&, cl::Buffer&,
	    cl::Kernel&, cl::Kernel&,
	    cl::CommandQueue&, cl::Context,
	    unsigned int*);

	void simulate_single_frame(
		particle*, particle*,
	    cl::Buffer&, cl::Buffer&,
	    unsigned int*,
        cl::Kernel&, cl::Kernel&, cl::Kernel&, cl::Kernel&,
	    cl::Kernel&, cl::CommandQueue&, cl::Context);
};

#endif