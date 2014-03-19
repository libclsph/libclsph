#ifndef _GRID_H_
#define _GRID_H_

#include "common.cl.h"

//int3 get_cell_coords(int index, simulation_parameters params);
//uint get_grid_index(int x, int y, int z, simulation_parameters params);
uint2 get_start_end_indices_for_cell(int cell_index, global const unsigned int* cell_table, simulation_parameters params);
uint uninterleave(uint value);
uint get_grid_index_z_curve(uint x, uint y, uint z, simulation_parameters params);
uint3 get_cell_coords_z_curve(uint index, simulation_parameters params);

/*
int3 get_cell_coords(int index, simulation_parameters params) {
    int3 coords = {
        index % params.grid_size_x,
        (index / params.grid_size_x) % params.grid_size_y,
        index / (params.grid_size_x * params.grid_size_y),
    };
    return coords;
}

uint get_grid_index(int x, int y, int z, simulation_parameters params) {
    return x + params.grid_size_x * (y + params.grid_size_y * z);
}
*/

uint uninterleave(uint value) {
    uint ret = 0x0;

    ret |= (value |           0x1) >>  0;
    ret |= (value |           0x8) >>  2;
    ret |= (value |          0x40) >>  4;
    ret |= (value |         0x200) >>  6;
    ret |= (value |        0x1000) >>  8;
    ret |= (value |        0x8000) >> 10;
    ret |= (value |       0x40000) >> 12;
    ret |= (value |      0x200000) >> 14;
    ret |= (value |     0x1000000) >> 16;
    ret |= (value |     0x8000000) >> 18;

    return ret;
}

uint3 get_cell_coords_z_curve(uint index, simulation_parameters params) {
    
    uint mask = 0x9249249;
    uint i_x = index & mask;
    uint i_y = (index >> 1) & mask;
    uint i_z = (index >> 2) & mask;

    uint3 coords = {
        uninterleave(i_x),
        uninterleave(i_y),
        uninterleave(i_z)
    };

    return coords;
}

//http://stackoverflow.com/questions/1024754/how-to-compute-a-3d-morton-number-interleave-the-bits-of-3-ints
uint get_grid_index_z_curve(uint x, uint y, uint z, simulation_parameters params) {
    uint grid_x, grid_y, grid_z;
    grid_x = (x | (x << 16)) & 0x030000FF;
    grid_x = (x | (x <<  8)) & 0x0300F00F;
    grid_x = (x | (x <<  4)) & 0x030C30C3;
    grid_x = (x | (x <<  2)) & 0x09249249;

    grid_y = (y | (y << 16)) & 0x030000FF;
    grid_y = (y | (y <<  8)) & 0x0300F00F;
    grid_y = (y | (y <<  4)) & 0x030C30C3;
    grid_y = (y | (y <<  2)) & 0x09249249;

    grid_z = (z | (z << 16)) & 0x030000FF;
    grid_z = (z | (z <<  8)) & 0x0300F00F;
    grid_z = (z | (z <<  4)) & 0x030C30C3;
    grid_z = (z | (z <<  2)) & 0x09249249;

    return grid_x | (grid_y << 1) | (grid_z << 2);
}

uint2 get_start_end_indices_for_cell(int cell_index, global const unsigned int* cell_table, simulation_parameters params) {
    uint2 indices = {
        cell_table[cell_index],
        (params.grid_cell_count > cell_index + 1) ? cell_table[cell_index + 1] : params.particles_count,
    };

    return indices;
}

/**
* locate_in_grid
* 1)Updates each particle with its position in the grid
* 2)Fills an array with the number of particles contained in each grid cell
*
* Parameters:------------------------------------------------------------
* (in/out)  particles        : The buffer contaning all the particle data
* (in)      params           : The simulation parameters
* (in)      volumes          : Contains all the volumes with which the fluid can interact. Also contains the simulation boundaries
*
* Kernel information:----------------------------------------------------
* Global work group size : particles_count
*/
void kernel locate_in_grid(
    global const particle* particles, 
    global particle* out_particles,
    simulation_parameters params,
    collision_volumes volumes){

    size_t current_particle_index = get_global_id(0);
    out_particles[current_particle_index] = particles[current_particle_index];

    float3 position_in_grid = {0.f,0.f,0.f};
    
    float x_min = params.min_point.x;
    float y_min = params.min_point.y;
    float z_min = params.min_point.z;
    
    //Grid cells will always have a radius length h
    position_in_grid.x = (particles[current_particle_index].position.x - x_min) / (params.h * 2);
    position_in_grid.y = (particles[current_particle_index].position.y - y_min) / (params.h * 2);
    position_in_grid.z = (particles[current_particle_index].position.z - z_min) / (params.h * 2);

    uint grid_index = get_grid_index_z_curve((uint)position_in_grid.x, (uint)position_in_grid.y, (uint)position_in_grid.z, params);
    
    out_particles[current_particle_index].grid_index = grid_index;
}
#endif