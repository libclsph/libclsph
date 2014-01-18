#define EXIT_ON_CL_ERROR

#include <iostream>
#include <string>

#include "sph_simulation.h"
#include "file_save_delegates/houdini_file_saver.h"

int main(int argc, char** argv) {

    if(argc < 5) {
        std::cout << "Too few arguments" << std::endl <<
            "Usage: ./sph <fluid_name> <simulation_properties_name> <scene_name> <frames_folder_prefix>" << std::endl;
        return -1;
    }

    sph_simulation simulation;
    houdini_file_saver saver = houdini_file_saver(std::string(argv[4]));

    try{
        simulation.load_settings(
            std::string("fluid_properties/") + argv[1] + std::string(".json"),
            std::string("simulation_properties/") + argv[2] + std::string(".json"));
    }
    catch (const std::exception& ex) {
        std::cerr << ex.what() << std::endl;
        std::exit(-1);
    }

    int i = 0;
    simulation.post_frame = [&] (particle* particles, const simulation_parameters& params, bool full_frame) {
        if(simulation.write_intermediate_frames != full_frame) {
            saver.writeFrameToFile(particles, params);
        }
        std::cout << "frame: " << i++ << std::endl;
    };

    std::cout << std::endl << 
        "Loaded parameters          " << std::endl <<
        "-----------------          " << std::endl <<
        "Simulation time:           " << simulation.parameters.simulation_time << std::endl <<
        "Target FPS:                " << simulation.parameters.target_fps << std::endl <<
        "Time delta:                " << simulation.parameters.time_delta << std::endl <<
        "Simulation scale:          " << simulation.parameters.simulation_scale << std::endl <<
        "Write intermediate frames: " << (simulation.write_intermediate_frames ? "true" : "false") << std::endl <<
        std::endl <<
        "Particle count:            " << simulation.parameters.particles_count << std::endl <<
        "Particle mass:             " << simulation.parameters.particle_mass << std::endl <<
        "Total mass:                " << simulation.parameters.total_mass << std::endl <<
        "Initial volume:            " << simulation.parameters.initial_volume << std::endl <<
        std::endl <<
        "Fluid density:             " << simulation.parameters.fluid_density << std::endl <<
        "Dynamic viscosity:         " << simulation.parameters.dynamic_viscosity << std::endl <<
        "Surface tension threshold: " << simulation.parameters.surface_tension_threshold << std::endl <<
        "Surface tension:           " << simulation.parameters.surface_tension << std::endl <<
        "Stiffness (k):             " << simulation.parameters.K << std::endl <<
        "Restitution:               " << simulation.parameters.restitution << std::endl <<
        std::endl <<
        "Kernel support radius (h): " << simulation.parameters.h << std::endl <<
        std::endl <<
        "Saving to folder:          " << saver.frames_folder_prefix + "frames/" << std::endl;

    simulation.load_scene(std::string("scenes/") + argv[3] + std::string(".json"));

    std::cout << std::endl << 
        "Revise simulation parameters.  Press q to quit, any other key to proceed with simulation" << std::endl;

    char response;
    std::cin >> response;
    if(response != 'q') {
        simulation.simulate();
    }

    return 0;
}
