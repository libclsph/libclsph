#define _USE_MATH_DEFINES
#include <cmath>

#define EXIT_ON_CL_ERROR

#include "../util/pico_json/picojson.h"
#include "../util/cereal/archives/binary.hpp"

#include "sph_simulation.h"
#include "common/util.h"

#define SORT_THREAD_COUNT 128
#define BUCKET_COUNT 256
#define RADIX_WIDTH 8

#define CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE_AMD 64

#define STEP_1 1
#define STEP_2 2
#define SORT_GPU 3
#define MEMORY_TRANSFERS 4
#define FILE_IO 5
#define SORT_CPU 6

const std::string BUFFER_KERNEL_FILE_NAME = "kernels/sph.cl";

cl::Device* running_device;

cl::Kernel make_kernel(cl::Program& p, const char* name) {

    cl_int cl_error; 
    cl::Kernel k = cl::Kernel(p, name, &cl_error);
    check_cl_error(cl_error);

    return k;
}

void set_kernel_args_internal(int, cl::Kernel&) {
    return;
}

template<typename T1, typename... TArgs>
void set_kernel_args_internal(int index, cl::Kernel& kernel, T1 a, TArgs... args) {
    check_cl_error(kernel.setArg(index, a));
    set_kernel_args_internal(index + 1, kernel, args...);
}

template<typename... TArgs>
void set_kernel_args(cl::Kernel& kernel, TArgs... args) {
    set_kernel_args_internal(0, kernel, args...);
}

/**
 * @brief Restores the state of the last simulation or places the particles in the shape of a cube if no state is found.
 *
 * @param[out] buffer        Points to the particles buffer to be filled
 * @param[in] parameters    Contains the simulation parameters
 *
 */
void sph_simulation::init_particles(particle* buffer , const simulation_parameters &parameters ) {

    int particles_per_cube_side = ceil(cbrtf(parameters.particles_count));
    float side_length = cbrtf(initial_volume);
    float spacing = side_length / (float)particles_per_cube_side;

    std::cout << "volume: " << initial_volume <<
    " side_length: " << side_length << " spacing: " << spacing << std::endl;

    //Last simualtion serialized its last frame
    //Lets load that and pick up where it let off
    std::filebuf fb;
    if (fb.open ("last_frame.bin",std::ios::in) )
    {
        std::istream file_in(&fb);

        cereal::BinaryInputArchive archive(file_in);
        archive.loadBinary(buffer,sizeof(particle)*parameters.particles_count);

        fb.close();
    }
    //No serialized particle array was found
    //Initialize the particles in their default position
    else{
        for(unsigned int i = 0; i < parameters.particles_count; ++i) {
            //Arrange the particles in the form of a cube
            buffer[i].position.s[0] = (float)(i % particles_per_cube_side) * spacing - side_length / 2.f;
            buffer[i].position.s[1] = ((float)((i / particles_per_cube_side) % particles_per_cube_side) * spacing);
            buffer[i].position.s[2] = (float)(i / (particles_per_cube_side * particles_per_cube_side)) * spacing - side_length / 2.f;

            buffer[i].velocity.s[0] = 0.f;
            buffer[i].velocity.s[1] = 0.f;
            buffer[i].velocity.s[2] = 0.f;
            buffer[i].intermediate_velocity.s[0] = 0.f;
            buffer[i].intermediate_velocity.s[1] = 0.f;
            buffer[i].intermediate_velocity.s[2] = 0.f;
            buffer[i].constant_acceleration.s[0] = parameters.constant_acceleration.s[0];
            buffer[i].constant_acceleration.s[1] = parameters.constant_acceleration.s[1];
            buffer[i].constant_acceleration.s[2] = parameters.constant_acceleration.s[2];

            buffer[i].density = 0.f;
            buffer[i].pressure = 0.f;
        }
    }

}

void sph_simulation::sort_particles(
    particle* particles,
    cl::Buffer& first_buffer, cl::Buffer& second_buffer,
    cl::Kernel& kernel_sort_count, cl::Kernel& kernel_sort,
    cl::CommandQueue& queue, cl::Context context,
    unsigned int* cell_table) {

    unsigned int* count_array = new unsigned int[SORT_THREAD_COUNT * BUCKET_COUNT];

    cl::Buffer count_buffer(
        context,  
        CL_MEM_READ_WRITE |  CL_MEM_ALLOC_HOST_PTR ,  
        sizeof(unsigned int) * SORT_THREAD_COUNT * BUCKET_COUNT);

    cl::Buffer* current_input_buffer = &first_buffer;
    cl::Buffer* current_output_buffer = &second_buffer;

    profile(MEMORY_TRANSFERS, &profile) {
        check_cl_error(
            queue.enqueueWriteBuffer(*current_input_buffer, CL_TRUE, 0,
                sizeof(particle) * parameters.particles_count, particles));
    }

    for(int pass_number = 0; pass_number < 4; ++pass_number) {
        for(int i = 0; i < SORT_THREAD_COUNT * BUCKET_COUNT; ++i) {
            count_array[i] = 0;
        }

        profile(MEMORY_TRANSFERS, &profile) {
            check_cl_error(
                queue.enqueueWriteBuffer(count_buffer, CL_TRUE, 0,
                    sizeof(unsigned int) * SORT_THREAD_COUNT * BUCKET_COUNT, count_array));
            queue.finish();
        }

        set_kernel_args(kernel_sort_count, 
            *current_input_buffer, count_buffer, parameters, SORT_THREAD_COUNT, pass_number, RADIX_WIDTH);

        profile(SORT_GPU, &profile) {
            check_cl_error(
                queue.enqueueNDRangeKernel(
                    kernel_sort_count, cl::NullRange,
                    cl::NDRange(SORT_THREAD_COUNT), cl::NullRange));
            queue.finish();
        }

        profile(MEMORY_TRANSFERS, &profile) {
            check_cl_error(
                queue.enqueueReadBuffer(
                    count_buffer, CL_TRUE, 0,
                    sizeof(unsigned int) * SORT_THREAD_COUNT * BUCKET_COUNT,
                    count_array));
            queue.finish();
        }

        profile(SORT_CPU, &profile) {
            unsigned int running_count = 0;
            for(int i = 0; i < SORT_THREAD_COUNT * BUCKET_COUNT; ++i) {
                unsigned int tmp = count_array[i];
                count_array[i] = running_count;
                running_count += tmp;
            }
        }

        profile(MEMORY_TRANSFERS, &profile) {
            check_cl_error(
                queue.enqueueWriteBuffer(count_buffer, CL_TRUE, 0,
                    sizeof(unsigned int) * SORT_THREAD_COUNT * BUCKET_COUNT,
                    count_array));
            queue.finish();
        }

        set_kernel_args(kernel_sort, 
            *current_input_buffer, *current_output_buffer, count_buffer, parameters, SORT_THREAD_COUNT, pass_number, RADIX_WIDTH);

        profile(SORT_GPU, &profile) {
            check_cl_error(
                queue.enqueueNDRangeKernel(
                    kernel_sort, cl::NullRange,
                    cl::NDRange(SORT_THREAD_COUNT), cl::NullRange));
            queue.finish();
        }

        cl::Buffer* tmp = current_input_buffer;
        current_input_buffer = current_output_buffer;
        current_output_buffer = tmp;
    }

    profile(MEMORY_TRANSFERS, &profile) {
        check_cl_error(
            queue.enqueueReadBuffer(*current_input_buffer, CL_TRUE, 0,
                sizeof(particle) * parameters.particles_count,
                particles));
        queue.finish();
    }

    delete[] count_array;

    profile(SORT_CPU, &profile) {
        unsigned int current_index = 0;
        for(unsigned int i = 0; i < parameters.grid_cell_count; ++i) {
            cell_table[i] = current_index;
            while(current_index != parameters.particles_count && particles[current_index].grid_index == i) {
                current_index++;
            }
        }
    }
}

void sph_simulation::simulate_single_frame(
    particle* in_particles, particle* out_particles,
    cl::Buffer& input_buffer, cl::Buffer& output_buffer,
    cl::Kernel& kernel_locate_in_grid, cl::Kernel& kernel_step_1,
    cl::Kernel& kernel_step_2, 
    cl::Kernel& kernel_sort_count, cl::Kernel& kernel_sort,
    cl::CommandQueue& queue, cl::Context context) {

    //-----------------------------------------------------
    // Calculate the optimal size for workgroups
    //-----------------------------------------------------

    //Start groups size at their maximum, make them smaller if necessary
    //Optimally parameters.particles_count should be devisible by CL_DEVICE_MAX_WORK_GROUP_SIZE
    //Refer to CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE

    unsigned int size_of_groups = running_device->getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>();
    while(parameters.particles_count%size_of_groups!=0){
        size_of_groups /= 2;
    }

    //Make sure that the workgroups are small enough and that the particle data will fit in local memory
    assert( size_of_groups <= running_device->getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>() );
    assert( size_of_groups*sizeof(particle) <= running_device->getInfo<CL_DEVICE_LOCAL_MEM_SIZE>() );

    //-----------------------------------------------------
    // Initial transfer to the GPU
    //-----------------------------------------------------

    profile(MEMORY_TRANSFERS, &profile) {
        check_cl_error(
            queue.enqueueWriteBuffer(
                input_buffer, CL_TRUE, 0,
                sizeof(particle) * parameters.particles_count, in_particles));
            queue.finish();
    }

    cl_float min_x, max_x, min_y, max_y, min_z, max_z;
    min_x = min_y = min_z = std::numeric_limits<cl_int>::max();
    max_x = max_y = max_z = std::numeric_limits<cl_int>::min();
    for(size_t i = 0; i < parameters.particles_count; ++i) {
        cl_float3 pos = in_particles[i].position;

        if(pos.s[0] < min_x) min_x = pos.s[0];
        if(pos.s[1] < min_y) min_y = pos.s[1];
        if(pos.s[2] < min_z) min_z = pos.s[2];

        if(pos.s[0] > max_x) max_x = pos.s[0];
        if(pos.s[1] > max_y) max_y = pos.s[1];
        if(pos.s[2] > max_z) max_z = pos.s[2];
    }

    parameters.min_point.s[0] = min_x;
    parameters.min_point.s[1] = min_y;
    parameters.min_point.s[2] = min_z;
    parameters.max_point.s[0] = max_x; 
    parameters.max_point.s[1] = max_y; 
    parameters.max_point.s[2] = max_z;

    parameters.grid_size_x = (cl_uint)((max_x - min_x) / (parameters.h * 2)) + 1;
    parameters.grid_size_y = (cl_uint)((max_y - min_y) / (parameters.h * 2)) + 1;
    parameters.grid_size_z = (cl_uint)((max_z - min_z) / (parameters.h * 2)) + 1;
    parameters.grid_cell_count = get_grid_index_z_curve(
        parameters.grid_size_x, 
        parameters.grid_size_y, 
        parameters.grid_size_z);

    unsigned int* cell_table = new unsigned int[parameters.grid_cell_count];

    //----------------------------------------------------------------
    // Locate each particle in the grid and build the grid count table
    //----------------------------------------------------------------

    set_kernel_args(kernel_locate_in_grid, input_buffer, output_buffer, parameters);

    check_cl_error(
        queue.enqueueNDRangeKernel(
            kernel_locate_in_grid, cl::NullRange,
            cl::NDRange(parameters.particles_count), cl::NDRange(size_of_groups)));

    profile(MEMORY_TRANSFERS, &profile) {
        check_cl_error(
            queue.enqueueReadBuffer(
                output_buffer,
                CL_TRUE, 0,
                sizeof(particle) * parameters.particles_count,
                out_particles));
            queue.finish();
    }

    sort_particles(
        out_particles,
        input_buffer, output_buffer,
        kernel_sort_count, kernel_sort,
        queue, context,
        cell_table);

    cl::Buffer cell_table_buffer(
        context,  
        CL_MEM_READ_WRITE |  CL_MEM_ALLOC_HOST_PTR ,  
        sizeof(unsigned int) * parameters.grid_cell_count);
    
    profile(MEMORY_TRANSFERS, &profile) {
        check_cl_error(
            queue.enqueueWriteBuffer(
                cell_table_buffer, CL_TRUE, 0,
                sizeof(unsigned int) * parameters.grid_cell_count, cell_table));

        check_cl_error(
            queue.enqueueWriteBuffer(
                input_buffer, CL_TRUE, 0,
                sizeof(particle) * parameters.particles_count, out_particles));
        queue.finish();
    }

    //-----------------------------------------------------
    // step_1
    //-----------------------------------------------------
    check_cl_error(kernel_step_1.setArg(0, input_buffer));
    check_cl_error(kernel_step_1.setArg(1, size_of_groups * sizeof(particle), NULL)); //Declare local memory in arguments
    check_cl_error(kernel_step_1.setArg(2, output_buffer));
    check_cl_error(kernel_step_1.setArg(3, parameters));
    check_cl_error(kernel_step_1.setArg(4, precomputed_terms));
    check_cl_error(kernel_step_1.setArg(5, cell_table_buffer));

    profile(STEP_1, &profile) {
        check_cl_error(
            queue.enqueueNDRangeKernel(
                kernel_step_1, cl::NullRange,
                cl::NDRange(parameters.particles_count), cl::NDRange(size_of_groups)));
                queue.finish();
    }

    profile(MEMORY_TRANSFERS, &profile) {
        check_cl_error(
            queue.enqueueReadBuffer(
                output_buffer,
                CL_TRUE, 0,
                sizeof(particle) * parameters.particles_count,
                out_particles));
        queue.finish();
    }

    //-----------------------------------------------------
    // step_2
    //-----------------------------------------------------
    cl::Buffer face_normals_buffer(context, CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR, sizeof(float) * current_scene.face_normals.size());
    cl::Buffer vertices_buffer(context, CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR, sizeof(float) * current_scene.vertices.size());
    cl::Buffer indices_buffer(context, CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR, sizeof(unsigned int) * current_scene.indices.size());

    profile(MEMORY_TRANSFERS, &profile) {
            check_cl_error(
                queue.enqueueWriteBuffer(
                    face_normals_buffer, CL_TRUE, 0,
                    sizeof(float) * current_scene.face_normals.size(), current_scene.face_normals.data()));

            check_cl_error(
                queue.enqueueWriteBuffer(
                    vertices_buffer, CL_TRUE, 0,
                    sizeof(float) * current_scene.vertices.size(), current_scene.vertices.data()));

            check_cl_error(
                queue.enqueueWriteBuffer(
                    indices_buffer, CL_TRUE, 0,
                    sizeof(unsigned int) * current_scene.indices.size(), current_scene.indices.data()));

            queue.finish();
    }

    check_cl_error(kernel_step_2.setArg(0, output_buffer));
    check_cl_error(kernel_step_2.setArg(1, input_buffer));
    check_cl_error(kernel_step_2.setArg(2, parameters));
    check_cl_error(kernel_step_2.setArg(3, precomputed_terms));
    check_cl_error(kernel_step_2.setArg(4, cell_table_buffer));
    check_cl_error(kernel_step_2.setArg(5, face_normals_buffer));
    check_cl_error(kernel_step_2.setArg(6, vertices_buffer));
    check_cl_error(kernel_step_2.setArg(7, indices_buffer));
    check_cl_error(kernel_step_2.setArg(8, current_scene.face_count));

    profile(STEP_2, &profile) {
        check_cl_error(
            queue.enqueueNDRangeKernel(
                kernel_step_2, cl::NullRange,
                cl::NDRange(parameters.particles_count), cl::NDRange(size_of_groups)));
        queue.finish();
    }

    profile(MEMORY_TRANSFERS, &profile) {
        check_cl_error(
            queue.enqueueReadBuffer(
                input_buffer, CL_TRUE, 0,
                sizeof(particle) * parameters.particles_count, out_particles));
        queue.finish();
    }

    delete[] cell_table;
}

void sph_simulation::simulate(int frame_count) {
    if(frame_count == 0) {
        frame_count = (int)ceil(parameters.simulation_time * parameters.target_fps);
    }

    cl_int cl_error; 

    cl::Context context;
    std::vector<cl::Device> device_array;
    check_cl_error(init_cl_single_device(&context, device_array, "", "", true));

    cl::CommandQueue queue(context, device_array[0], 0, &cl_error);
    check_cl_error(cl_error);

    running_device = &device_array[0];

    std::string source = readKernelFile(BUFFER_KERNEL_FILE_NAME);
    cl::Program program;
    check_cl_error(make_program(&program, context, device_array, source, true, "-I ./kernels/ -I ./common/"));

    cl::Kernel kernel_step_1 = make_kernel(program, "step_1");
    cl::Kernel kernel_step_2 = make_kernel(program, "step_2");
    cl::Kernel kernel_locate_in_grid = make_kernel(program, "locate_in_grid");
    cl::Kernel kernel_sort_count = make_kernel(program, "sort_count");
    cl::Kernel kernel_sort = make_kernel(program, "sort");

    cl::Buffer input_buffer(context,  CL_MEM_READ_WRITE |  CL_MEM_ALLOC_HOST_PTR ,  sizeof(particle) * parameters.particles_count);
    cl::Buffer output_buffer(context, CL_MEM_READ_WRITE |  CL_MEM_ALLOC_HOST_PTR ,  sizeof(particle) * parameters.particles_count);

    particle* particles = new particle[parameters.particles_count];
    init_particles(particles,parameters);

    std::cout << std::endl;

    for(int i = 0; i < frame_count; ++i) {
        if(pre_frame)  {
            pre_frame(particles, parameters, true, profile);
        }

        for(int j = 0; (float)j < (1.f / parameters.simulation_scale); ++j) {
            if(pre_frame) pre_frame(particles, parameters, false, profile);

            simulate_single_frame(
                particles, particles,
                input_buffer, output_buffer,
                kernel_locate_in_grid, kernel_step_1, kernel_step_2,
                kernel_sort_count, kernel_sort,
                queue, context);

            if(post_frame) post_frame(particles, parameters, false, profile);
        }

        if(post_frame) {
            post_frame(particles, parameters, true, profile);
        }
    }
    delete[] particles;
}

void sph_simulation::load_settings(std::string fluid_file_name, std::string parameters_file_name) {
    int particles_inside_influence_radius = 0;

    {
        picojson::value fluid_params;
        std::ifstream fluid_stream(fluid_file_name);
        fluid_stream >> fluid_params;

        parameters.fluid_density = (float)(fluid_params.get<picojson::object>()["fluid_density"].get<double>());
        parameters.dynamic_viscosity = (float)(fluid_params.get<picojson::object>()["dynamic_viscosity"].get<double>());
        parameters.restitution = (float)(fluid_params.get<picojson::object>()["restitution"].get<double>());
        if (parameters.restitution < 0 || parameters.restitution > 1) {
            throw std::runtime_error( "Restitution has an invalid value!" );
        }

        parameters.K = (float)(fluid_params.get<picojson::object>()["k"].get<double>());
        parameters.surface_tension_threshold = (float)(fluid_params.get<picojson::object>()["surface_tension_threshold"].get<double>());
        parameters.surface_tension = (float)(fluid_params.get<picojson::object>()["surface_tension"].get<double>());
        particles_inside_influence_radius = (int)(fluid_params.get<picojson::object>()["particles_inside_influence_radius"].get<double>());
    }

    {
        picojson::value sim_params;
        std::ifstream sim_stream(parameters_file_name);
        sim_stream >> sim_params;

        parameters.particles_count = (unsigned int)(sim_params.get<picojson::object>()["particles_count"].get<double>());

        if( parameters.particles_count%CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE_AMD != 0 ){
            std::cout << std::endl << "\033[1;31m You should choose a number of particles that is divisble by the preferred work group size.\033[0m";
            std::cout << std::endl << "\033[1;31m Performances will be sub-optimal.\033[0m" << std::endl;
        }

        parameters.particle_mass = (float)(sim_params.get<picojson::object>()["particle_mass"].get<double>());
        parameters.simulation_time = (float)(sim_params.get<picojson::object>()["simulation_time"].get<double>());
        parameters.target_fps = (float)(sim_params.get<picojson::object>()["target_fps"].get<double>());
        parameters.simulation_scale = (float)(sim_params.get<picojson::object>()["simulation_scale"].get<double>());

        parameters.constant_acceleration.s[0] = (float)(sim_params.get<picojson::object>()["constant_acceleration"].get<picojson::object>()["x"].get<double>());
        parameters.constant_acceleration.s[1] = (float)(sim_params.get<picojson::object>()["constant_acceleration"].get<picojson::object>()["y"].get<double>());
        parameters.constant_acceleration.s[2] = (float)(sim_params.get<picojson::object>()["constant_acceleration"].get<picojson::object>()["z"].get<double>());

        write_intermediate_frames = sim_params.get<picojson::object>()["write_all_frames"].get<bool>();
        serialize = sim_params.get<picojson::object>()["serialize"].get<bool>();
    }

    parameters.total_mass = parameters.particles_count * parameters.particle_mass;
    initial_volume = parameters.total_mass / parameters.fluid_density;
    parameters.h = cbrtf(
        3.f * 
        (particles_inside_influence_radius * 
            (initial_volume / parameters.particles_count)) / 
        (4.f * M_PI));
    parameters.time_delta = 1.f / parameters.target_fps;

    parameters.max_velocity = 0.8f * parameters.h / parameters.time_delta;

    precomputed_terms.poly_6 = 315.f / (64.f * M_PI * pow(parameters.h, 9.f));
    precomputed_terms.poly_6_gradient = -945.f / (32.f * M_PI * pow(parameters.h, 9.f));
    precomputed_terms.poly_6_laplacian = -945.f / (32.f * M_PI * pow(parameters.h, 9.f));
    precomputed_terms.spiky = -45.f / (M_PI * pow(parameters.h, 6.f));
    precomputed_terms.viscosity = 45.f / (M_PI * pow(parameters.h, 6.f));
}
