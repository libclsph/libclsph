#define EXIT_ON_CL_ERROR

#include <iostream>
#include <iomanip>
#include <string>

#include "sph_simulation.h"
#include "file_save_delegates/houdini_file_saver.h"
#include "util/cereal/archives/binary.hpp"
#include "util/profile.hpp"

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

    simulation.pre_frame = [&] (particle* particles, const simulation_parameters& params, bool full_frame, profile_data& profile) {
        if(simulation.write_intermediate_frames != full_frame) {
            saver.writeFrameToFile(particles, params);
        }

        //Lets serialize the frames if it was specified
        if(simulation.serialize){
            std::ofstream file_out( "last_frame.bin" );
            cereal::BinaryOutputArchive archive(file_out);
            archive.saveBinary(particles,sizeof(particle)*params.particles_count);
        }

        ++i;

        std::cout << "Profile values averages over all frames" << std::endl;

        std::cout << std::setiosflags(std::ios::fixed) <<
        "|" << "            step 1"
        "|" << "            step 2"
        "|" << "   sort (GPU part)"  
        "|" << "   sort (CPU part)"  
        "|" << "  memory transfers"
        "|" << "           file IO" << 
        "|" << std::endl;

        std::cout << std::setiosflags(std::ios::fixed) <<
        "|" << std::setprecision(3) << std::setw(18) << (float)profile.data[1] / i / 1000.f <<
        "|" << std::setprecision(3) << std::setw(18) << (float)profile.data[2] / i / 1000.f <<
        "|" << std::setprecision(3) << std::setw(18) << (float)profile.data[3] / i / 1000.f <<
        "|" << std::setprecision(3) << std::setw(18) << (float)profile.data[6] / i / 1000.f <<
        "|" << std::setprecision(3) << std::setw(18) << (float)profile.data[4] / i / 1000.f <<
        "|" << std::setprecision(3) << std::setw(18) << (float)profile.data[5] / i / 1000.f << 
        "|" << std::endl;

        int num_frames = 
            params.simulation_time / (params.time_delta * params.simulation_scale);

        int progress = i * 80 / num_frames;

        std::cout << std::endl << "Progress" << std::endl;
        std::cout << "[";
        for(int j = 0; j < 80; ++j) {
            if(j < progress) {
                std::cout << "-";
            } else {
                std::cout << " ";
            }
        }
        std::cout << "] " << i << "/" << num_frames << std::endl;
        std::cout << "\033[F" << "\033[F" << "\033[F" << "\033[F" << "\033[F" << "\033[F";
    };

    std::cout << std::endl << 
        "Loaded parameters          " << std::endl <<
        "-----------------          " << std::endl <<
        "Simulation time:           " << simulation.parameters.simulation_time << std::endl <<
        "Target FPS:                " << simulation.parameters.target_fps << std::endl <<
        "Time delta:                " << simulation.parameters.time_delta << std::endl <<
        "Simulation scale:          " << simulation.parameters.simulation_scale << std::endl <<
        "Write intermediate frames: " << (simulation.write_intermediate_frames ? "true" : "false") << std::endl <<
        "Serialize frames:          " << (simulation.serialize ? "true" : "false") << std::endl <<
        std::endl <<
        "Particle count:            " << simulation.parameters.particles_count << std::endl <<
        "Particle mass:             " << simulation.parameters.particle_mass << std::endl <<
        "Total mass:                " << simulation.parameters.total_mass << std::endl <<
        "Initial volume:            " << simulation.initial_volume << std::endl <<
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

    if(!simulation.current_scene.load(argv[3])) {
        std::cerr << "Unable to load scene: " << argv[3] << std::endl;
        return -1;
    }

    //If the serialization data is not the right size, delete it
    //This probably means the last simulation ran with a different number of particles or the serialization was interrupted
    std::filebuf fb;
    if (fb.open ("last_frame.bin",std::ios::in))
    {
        std::istream file_in(&fb);

        file_in.seekg(0,std::ios_base::end);

        size_t file_size = file_in.tellg();

        //This information is important, it indicates that the behavior of the simulator is completely different, use color to draw attention of user

        if(file_size == simulation.parameters.particles_count*sizeof(particle)){
            std::cout << std::endl << "\033[1;32m Serialized frame found. " << " Simulation will pick up where last run left off.\033[0m";
            std::cout << std::endl << "\033[1;32m To start a new simulation, delete last_frame.bin. \033[0m" << std::endl;
        }
        else{
            std::cout << std::endl << "\033[1;31m Serialized frame of incorrect size found. Revert to last know settings or delete it, then try again. \033[0m" << std::endl;
            return 0;
        }

        fb.close();
    }

    std::cout << std::endl << 
        "Revise simulation parameters.  Press q to quit, any other key to proceed with simulation" << std::endl;

    char response;
    std::cin >> response;
    if(response != 'q') {
        simulation.simulate();
    }

    return 0;
}
