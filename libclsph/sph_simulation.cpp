#define _USE_MATH_DEFINES
#include <cmath>

#define EXIT_ON_CL_ERROR

#include "../util/pico_json/picojson.h"
#include "../util/profile.hpp"

#include "sph_simulation.h"
#include "collision_volumes_loader.h"
#include <assert.h>

const std::string BUFFER_KERNEL_FILE_NAME = "kernels/sph.cl";

cl::Device* running_device;

void sph_simulation::init_particles(particle* buffer, int count) {
    int particles_per_cube_side = ceil(cbrtf(parameters.particles_count));
    float side_length = cbrtf(parameters.initial_volume);
    float spacing = side_length / (float)particles_per_cube_side;

    std::cout << "volume: " << parameters.initial_volume <<
    " side_length: " << side_length << " spacing: " << spacing << std::endl;

    for(int i = 0; i < count; ++i) {
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

#define SORT_THREAD_COUNT 1024
#define BUCKET_COUNT 256
#define RADIX_WIDTH 8

void sph_simulation::sort_particles(
    particle* particles,
    cl::Buffer& input_buffer, cl::Buffer& output_buffer,
    cl::Kernel& kernel_sort_count, cl::Kernel& kernel_sort,
    cl::CommandQueue& queue, cl::Context context,
    unsigned int* cell_table) {

    particle* debug = new particle[parameters.particles_count];
    for(int i = 0; i < parameters.particles_count; ++i) {
        debug[i] = particles[i];
    }

    unsigned int* count_array = new unsigned int[SORT_THREAD_COUNT * BUCKET_COUNT];

    cl::Buffer count_buffer(
        context,  
        CL_MEM_READ_WRITE |  CL_MEM_ALLOC_HOST_PTR ,  
        sizeof(unsigned int) * SORT_THREAD_COUNT * BUCKET_COUNT);

    for(int pass_number = 0; pass_number < 4; ++pass_number) {

        for(int i = 0; i < SORT_THREAD_COUNT * BUCKET_COUNT; ++i) {
            count_array[i] = 0;
        }

        check_cl_error(
            queue.enqueueWriteBuffer(input_buffer, CL_TRUE, 0,
                sizeof(particle) * parameters.particles_count, particles));

        check_cl_error(
            queue.enqueueWriteBuffer(count_buffer, CL_TRUE, 0,
                sizeof(unsigned int) * SORT_THREAD_COUNT * BUCKET_COUNT, count_array));

        check_cl_error(kernel_sort_count.setArg(0, input_buffer));
        check_cl_error(kernel_sort_count.setArg(1, count_buffer));
        check_cl_error(kernel_sort_count.setArg(2, parameters));
        check_cl_error(kernel_sort_count.setArg(3, SORT_THREAD_COUNT));
        check_cl_error(kernel_sort_count.setArg(4, pass_number));
        check_cl_error(kernel_sort_count.setArg(5, RADIX_WIDTH));

        check_cl_error(
            queue.enqueueNDRangeKernel(
                kernel_sort_count, cl::NullRange,
                cl::NDRange(SORT_THREAD_COUNT), cl::NullRange));

        check_cl_error(
            queue.enqueueReadBuffer(
                count_buffer, CL_TRUE, 0,
                sizeof(unsigned int) * SORT_THREAD_COUNT * BUCKET_COUNT,
                count_array));

        unsigned int running_count = 0;
        for(int i = 0; i < SORT_THREAD_COUNT * BUCKET_COUNT; ++i) {
            unsigned int tmp = count_array[i];
            count_array[i] = running_count;
            running_count += tmp;
        }

        check_cl_error(
            queue.enqueueWriteBuffer(count_buffer, CL_TRUE, 0,
                sizeof(unsigned int) * SORT_THREAD_COUNT * BUCKET_COUNT,
                count_array));

        check_cl_error(kernel_sort.setArg(0, input_buffer));
        check_cl_error(kernel_sort.setArg(1, output_buffer));
        check_cl_error(kernel_sort.setArg(2, count_buffer));
        check_cl_error(kernel_sort.setArg(3, parameters));
        check_cl_error(kernel_sort.setArg(4, SORT_THREAD_COUNT));
        check_cl_error(kernel_sort.setArg(5, pass_number));
        check_cl_error(kernel_sort.setArg(6, RADIX_WIDTH));

        check_cl_error(
            queue.enqueueNDRangeKernel(
                kernel_sort, cl::NullRange,
                cl::NDRange(SORT_THREAD_COUNT), cl::NullRange));

        check_cl_error(
            queue.enqueueReadBuffer(output_buffer, CL_TRUE, 0,
                sizeof(particle) * parameters.particles_count,
                particles));
    }

    delete[] count_array;

    unsigned int current_index = 0;
    for(int i = 0; i < parameters.grid_cell_count; ++i) {
        cell_table[i] = current_index;
        while(current_index != parameters.particles_count && particles[current_index].grid_index == i) {
            current_index++;
        }
    }
}

void sph_simulation::simulate_single_frame(
    particle* in_particles, particle* out_particles,
    cl::Buffer& input_buffer, cl::Buffer& output_buffer,
    unsigned int* cell_table,
    cl::Kernel& kernel_locate_in_grid,cl::Kernel& kernel_build_grid, cl::Kernel& kernel_step_1,
    cl::Kernel& kernel_step_2, 
    cl::Kernel& kernel_sort_count, cl::Kernel& kernel_sort,
    cl::CommandQueue& queue, cl::Context context) {

    //-----------------------------------------------------
    // Calculate the optimal size for workgroups
    //-----------------------------------------------------

    //Start groups size at their maximum, make them smaller if necessary
    //Optimally parameters.particles_count should be devisible by CL_DEVICE_MAX_WORK_GROUP_SIZE
    //Refer to CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE
    //unsigned int work_group_size_multiple = kernel_locate_in_grid.getWorkGroupInfo<CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE>(*running_device);

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

    check_cl_error(
        queue.enqueueWriteBuffer(
            input_buffer, CL_TRUE, 0,
            sizeof(particle) * parameters.particles_count, in_particles));

    //----------------------------------------------------------------
    // Locate each particle in the grid and build the grid count table
    //----------------------------------------------------------------

    check_cl_error(kernel_locate_in_grid.setArg(0, input_buffer));
    check_cl_error(kernel_locate_in_grid.setArg(1, output_buffer));
    check_cl_error(kernel_locate_in_grid.setArg(2, parameters));
    check_cl_error(kernel_locate_in_grid.setArg(3, volumes));

    check_cl_error(
        queue.enqueueNDRangeKernel(
            kernel_locate_in_grid, cl::NullRange,
            cl::NDRange(parameters.particles_count), cl::NDRange(size_of_groups)));

    check_cl_error(
        queue.enqueueReadBuffer(
            output_buffer,
            CL_TRUE, 0,
            sizeof(particle) * parameters.particles_count,
            out_particles));

    sort_particles(
        out_particles,
        input_buffer, output_buffer,
        kernel_sort_count, kernel_sort,
        queue, context,
        cell_table);

    cl::Buffer cell_table_buffer(
        context,  
        CL_MEM_READ_WRITE |  CL_MEM_ALLOC_HOST_PTR ,  
        sizeof(unsigned int) * SORT_THREAD_COUNT * BUCKET_COUNT);
    
    check_cl_error(
        queue.enqueueWriteBuffer(
            cell_table_buffer, CL_TRUE, 0,
            sizeof(unsigned int) * parameters.grid_cell_count, cell_table));

    check_cl_error(
        queue.enqueueWriteBuffer(
            input_buffer, CL_TRUE, 0,
            sizeof(particle) * parameters.particles_count, out_particles));

    //-----------------------------------------------------
    // step_1
    //-----------------------------------------------------
            
    check_cl_error(kernel_step_1.setArg(0, input_buffer));
    check_cl_error(kernel_step_1.setArg(1, size_of_groups * sizeof(particle) , NULL)); //Declare local memory in arguments
    check_cl_error(kernel_step_1.setArg(2, output_buffer));
    check_cl_error(kernel_step_1.setArg(3, parameters));
    check_cl_error(kernel_step_1.setArg(4, cell_table_buffer));

    profile("execute range kernel (step 1)", profile_block::COUT_LOG) {

    check_cl_error(
        queue.enqueueNDRangeKernel(
            kernel_step_1, cl::NullRange, 
            cl::NDRange(parameters.particles_count), cl::NDRange(size_of_groups)));
            queue.finish();
    }

    check_cl_error(
        queue.enqueueReadBuffer(
            output_buffer,
            CL_TRUE, 0,
            sizeof(particle) * parameters.particles_count,
            out_particles));

    //-----------------------------------------------------
    // step_2
    //-----------------------------------------------------

    check_cl_error(
        queue.enqueueWriteBuffer(
            input_buffer, CL_TRUE, 0,
            sizeof(particle) * parameters.particles_count, out_particles));

    check_cl_error(kernel_step_2.setArg(0, input_buffer));
    check_cl_error(kernel_step_2.setArg(1, output_buffer));
    check_cl_error(kernel_step_2.setArg(2, parameters));
    check_cl_error(kernel_step_2.setArg(3, volumes));
    check_cl_error(kernel_step_2.setArg(4, cell_table_buffer));

    check_cl_error(
        queue.enqueueNDRangeKernel(
            kernel_step_2, cl::NullRange, 
            cl::NDRange(parameters.particles_count), cl::NDRange(size_of_groups)));

    //Read the result back from the device
    check_cl_error(
        queue.enqueueReadBuffer(
            output_buffer, CL_TRUE, 0,
            sizeof(particle) * parameters.particles_count, out_particles));
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

    cl::Kernel kernel_step_1 = cl::Kernel(program, "step_1", &cl_error);
    check_cl_error(cl_error);

    cl::Kernel kernel_step_2 = cl::Kernel(program, "step_2", &cl_error);
    check_cl_error(cl_error);

    cl::Kernel kernel_locate_in_grid = cl::Kernel(program, "locate_in_grid", &cl_error);
    check_cl_error(cl_error);

    cl::Kernel kernel_build_grid = cl::Kernel(program, "build_grid", &cl_error);
    check_cl_error(cl_error);

    cl::Kernel kernel_sort_count = cl::Kernel(program, "sort_count", &cl_error);
    check_cl_error(cl_error);

    cl::Kernel kernel_sort = cl::Kernel(program, "sort", &cl_error);
    check_cl_error(cl_error);

    parameters.grid_size_x = (int)(volumes.scene_bounding_box.axis_extends.s[0] * 2 / (parameters.h * 2));
    parameters.grid_size_y = (int)(volumes.scene_bounding_box.axis_extends.s[1] * 2 / (parameters.h * 2));
    parameters.grid_size_z = (int)(volumes.scene_bounding_box.axis_extends.s[2] * 2 / (parameters.h * 2));
    parameters.grid_cell_count = parameters.grid_size_x * parameters.grid_size_y * parameters.grid_size_z;

    std::cout << "Grid size:" << parameters.grid_size_x << ":" << parameters.grid_size_y << ":" << parameters.grid_size_z << std::endl;

    unsigned int* cell_table = new unsigned int[parameters.grid_cell_count];

    cl::Buffer input_buffer(context,  CL_MEM_READ_WRITE |  CL_MEM_ALLOC_HOST_PTR ,  sizeof(particle) * parameters.particles_count);
    cl::Buffer output_buffer(context, CL_MEM_READ_WRITE |  CL_MEM_ALLOC_HOST_PTR ,  sizeof(particle) * parameters.particles_count);

    particle* particles = new particle[parameters.particles_count];
    init_particles(particles, parameters.particles_count);

    std::cout << std::endl;

    for(int i = 0; i < frame_count; ++i) {
        if(pre_frame) pre_frame(particles, parameters, true);

        for(int j = 0; (float)j < (1.f / parameters.simulation_scale); ++j) {
            if(pre_frame) pre_frame(particles, parameters, false);

            simulate_single_frame(
                particles, particles,
                input_buffer, output_buffer,
                cell_table,
                kernel_locate_in_grid, kernel_build_grid, kernel_step_1, kernel_step_2, 
                kernel_sort_count, kernel_sort,
                queue, context);

            if(post_frame) post_frame(particles, parameters, false);
        }

        if(post_frame) post_frame(particles, parameters, true);
    }
    delete[] particles;
    delete[] cell_table;
}

void sph_simulation::load_scene(std::string scene_file_name) {
    picojson::value scene;
    std::ifstream scene_stream(scene_file_name);

    collision_volumes_loader loader;

    volumes = loader.load_standard_json(scene_stream);
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
        if ( parameters.restitution <= 0 || parameters.restitution > 1) {
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
        parameters.particle_mass = (float)(sim_params.get<picojson::object>()["particle_mass"].get<double>());
        parameters.simulation_time = (float)(sim_params.get<picojson::object>()["simulation_time"].get<double>());
        parameters.target_fps = (float)(sim_params.get<picojson::object>()["target_fps"].get<double>());
        parameters.simulation_scale = (float)(sim_params.get<picojson::object>()["simulation_scale"].get<double>());

        parameters.constant_acceleration.s[0] = (float)(sim_params.get<picojson::object>()["constant_acceleration"].get<picojson::object>()["x"].get<double>());
        parameters.constant_acceleration.s[1] = (float)(sim_params.get<picojson::object>()["constant_acceleration"].get<picojson::object>()["y"].get<double>());
        parameters.constant_acceleration.s[2] = (float)(sim_params.get<picojson::object>()["constant_acceleration"].get<picojson::object>()["z"].get<double>());

        write_intermediate_frames = sim_params.get<picojson::object>()["write_all_frames"].get<bool>();
    }

    parameters.total_mass = parameters.particles_count * parameters.particle_mass;
    parameters.initial_volume = parameters.total_mass / parameters.fluid_density;
    parameters.h = cbrtf(
        3.f * 
        (particles_inside_influence_radius * 
            (parameters.initial_volume / parameters.particles_count)) / 
        (4.f * M_PI));
    parameters.time_delta = 1.f / parameters.target_fps;
}
