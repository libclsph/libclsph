#ifndef _GRID_H_
#define _GRID_H_

#include "common.cl"

int3 get_cell_coords(int index, simulation_parameters params);
int get_grid_index(int x, int y, int z, simulation_parameters params);
uint2 get_start_end_indices_for_cell(int cell_index, global const unsigned int* cell_table, simulation_parameters params);

/**
 * @brief      Find the cartesian coordinates of a grid cell from its index in the flattened array
 *
 * @param[in]  index    The index of the grid cell in the flattened array
 * @param[in]  params   The simulation parameters
 *
 * @return     The cartesian coordinates of the grid cell_index
 *
 */
int3 get_cell_coords(int index, simulation_parameters params) {
    int3 coords = {
        index % params.grid_size_x,
        (index / params.grid_size_x) % params.grid_size_y,
        index / (params.grid_size_x * params.grid_size_y),
    };
    return coords;
}

/**
 * @brief      Finds the index of a grid cell in the flattened cell table from its cartesian coordinates
 *
 * @param[in]  x        The x position of the cell in the cartesian grid.
 * @param[in]  y        The y position of the cell in the cartesian grid.
 * @param[in]  z        The z position of the cell in the cartesian grid.
 * @param[in]  params   The simulation parameters
 *
 * @return     The index of the grid cell
 *
 */
int get_grid_index(int x, int y, int z, simulation_parameters params) {
    return x + params.grid_size_x * (y + params.grid_size_y * z);
}

/**
 * @brief      Locates the particle data for a certain grid cell in the cell table.
 *
 * @param[in]  cell_index       The index of the grid cell we are examining.
 * @param[in]  cell_table       A flattened representation of the grid contents
 * @param[in]  params           Contains the simulation parameters
 *
 * @return     The start and finish indexes of the subarray that contains the particles that can be found at cell_index.
 *
 */
uint2 get_start_end_indices_for_cell(int cell_index, global const unsigned int* cell_table, simulation_parameters params) {

    uint2 indices = {
        cell_table[cell_index],
        (params.grid_cell_count > cell_index + 1) ? cell_table[cell_index + 1] : params.particles_count,
    };

    return indices;
}

/**
 * @brief Updates each particle with its position in the grid and fills an array with the number of particles contained in each grid cell
 *
 * @param[in]  particles        Contains all the particle data
 * @param[out] out_particles    Will contain the particle data with the added information
 * @param[in]  params           Contains the simulation parameters
 * @param[in]  volumes          Contains all the volumes with which the fluid can interact. Also contains the simulation boundaries
 *
 */
void kernel locate_in_grid(
    global const particle* particles, 
    global particle* out_particles,
    simulation_parameters params,
    collision_volumes volumes){

    const size_t current_particle_index = get_global_id(0);
    out_particles[current_particle_index] = particles[current_particle_index];

    float3 position_in_grid = {0.f,0.f,0.f};
    
    float x_min = params.min_point.x;
    float y_min = params.min_point.y;
    float z_min = params.min_point.z;
    
    //Grid cells will always have a radius length h
    position_in_grid.x = (particles[current_particle_index].position.x - x_min) / (params.h * 2);
    position_in_grid.y = (particles[current_particle_index].position.y - y_min) / (params.h * 2);
    position_in_grid.z = (particles[current_particle_index].position.z - z_min) / (params.h * 2);

    int grid_index = get_grid_index((int)position_in_grid.x, (int)position_in_grid.y, (int)position_in_grid.z, params);
    
    out_particles[current_particle_index].grid_index = grid_index;
}
#endif
