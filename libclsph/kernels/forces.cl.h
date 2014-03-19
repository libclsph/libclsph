#include "common.cl.h"
#include "grid.h"

float compute_density(size_t current_particle_index, global const particle* others, const simulation_parameters params);
float3 compute_internal_forces(size_t current_particle_index, global const particle* others, const simulation_parameters params);
float compute_density_with_grid(size_t current_particle_index, global const particle* others, const simulation_parameters params, global const unsigned int* grid_cell_particle_list);
float3 compute_internal_forces_with_grid(size_t current_particle_index, global const particle* others, const simulation_parameters params, global const unsigned int* grid_cell_particle_list);

float compute_density_with_grid(
	size_t current_particle_index, 
	global const particle* others, 
	const simulation_parameters params,
	global const unsigned int* grid_cell_particle_list) {

	float density = 0.f;

	uint3 cell_coords = get_cell_coords_z_curve(others[current_particle_index].grid_index, params);

	for (int z = cell_coords.z - 1; z <= cell_coords.z + 1; ++z) {
    	if(z < 0 || z >= params.grid_size_z) continue;

        for (int y = cell_coords.y - 1; y <= cell_coords.y + 1; ++y) {
    		if(y < 0 || y >= params.grid_size_y) continue;

            for (int x = cell_coords.x - 1; x <= cell_coords.x + 1; ++x) {
    			if(x < 0 || x >= params.grid_size_x) continue;

    			uint grid_index = get_grid_index_z_curve(x, y, z, params);
    			uint2 indices = get_start_end_indices_for_cell(grid_index, grid_cell_particle_list, params);

    			for(size_t i = indices.x; i < indices.y; ++i) {

    				density += params.particle_mass * poly_6(
						distance(others[current_particle_index].position, others[i].position), 
						params.h);
    			}
    		}
    	}
    }

	return density;
}

float3 compute_internal_forces_with_grid(
	size_t current_particle_index, 
	global const particle* others, 
	const simulation_parameters params,
	global const unsigned int* grid_cell_particle_list) {

	float3 pressure_term = { 0.f, 0.f, 0.f };
	float3 viscosity_term = { 0.f, 0.f, 0.f };
	//compute the inward surface normal, it's the gradient of the color field
	float3 normal = { 0.f, 0.f, 0.f };
	//also need the color field laplacian
	float color_field_laplacian = 0.f;

	uint3 cell_coords = get_cell_coords_z_curve(others[current_particle_index].grid_index, params);

	for (int z = cell_coords.z - 1; z <= cell_coords.z + 1; ++z) {
    	if(z < 0 || z >= params.grid_size_z) continue;
        for (int y = cell_coords.y - 1; y <= cell_coords.y + 1; ++y) {
    		if(y < 0 || y >= params.grid_size_y) continue;
            for (int x = cell_coords.x - 1; x <= cell_coords.x + 1; ++x) {
    			if(x < 0 || x >= params.grid_size_x) continue;

    			uint grid_index = get_grid_index_z_curve(x, y, z, params);
    			uint2 indices = get_start_end_indices_for_cell(grid_index, grid_cell_particle_list, params);

    			for(size_t i = indices.x; i < indices.y; ++i) {

    				if(i != current_particle_index) {
						//[kelager] (4.11)
						pressure_term += 
							(others[i].pressure / pown(others[i].density, 2) + 
							others[current_particle_index].pressure / pown(others[current_particle_index].density, 2)) * 
							params.particle_mass * 
							spiky_gradient(others[current_particle_index].position - others[i].position, params.h);

						viscosity_term += 
							(others[i].velocity - others[current_particle_index].velocity) *
							(params.particle_mass / others[i].density) *
							viscosity_laplacian(
								length(others[current_particle_index].position - others[i].position),
								params.h);
					}

					normal += 
						params.particle_mass / others[i].density * 
						poly_6_gradient(others[current_particle_index].position - others[i].position, params.h);

					color_field_laplacian += 
						params.particle_mass / others[i].density * 
						poly_6_laplacian(length(others[current_particle_index].position - others[i].position), params.h);
    			}
    		}
    	}
    }

	float3 sum = 
		(-others[current_particle_index].density * pressure_term) + 
		(viscosity_term * params.dynamic_viscosity);

	if(length(normal) > params.surface_tension_threshold) {
		sum += -params.surface_tension * color_field_laplacian * normal / length(normal);
	}
	
	return sum;
}

float compute_density(size_t current_particle_index, global const particle* others, const simulation_parameters params) {
    float density = 0.f;

	for(size_t i = 0; i < params.particles_count; ++i) {
		density += params.particle_mass * poly_6(
			distance(others[current_particle_index].position, others[i].position), 
				params.h);
	}

	return density;
}

float3 compute_internal_forces(size_t current_particle_index, global const particle* others, const simulation_parameters params) {

	float3 pressure_term = { 0.f, 0.f, 0.f };
	float3 viscosity_term = { 0.f, 0.f, 0.f };
	//compute the inward surface normal, it's the gradient of the color field
	float3 normal = { 0.f, 0.f, 0.f };
	//also need the color field laplacian
	float color_field_laplacian = 0.f;

	for(size_t i = 0; i < params.particles_count; ++i) {
		if(i != current_particle_index) {
			//[kelager] (4.11)
			pressure_term += 
				(others[i].pressure / pown(others[i].density, 2) + 
				others[current_particle_index].pressure / pown(others[current_particle_index].density, 2)) * 
				params.particle_mass * 
				spiky_gradient(others[current_particle_index].position - others[i].position, params.h);

			viscosity_term += 
				(others[i].velocity - others[current_particle_index].velocity) *
				(params.particle_mass / others[i].density) *
				viscosity_laplacian(
					length(others[current_particle_index].position - others[i].position),
					params.h);
		}

		normal += 
			params.particle_mass / others[i].density * 
			poly_6_gradient(others[current_particle_index].position - others[i].position, params.h);

		color_field_laplacian += 
			params.particle_mass / others[i].density * 
			poly_6_laplacian(length(others[current_particle_index].position - others[i].position), params.h);
	}

	float3 sum = 
		(-others[current_particle_index].density * pressure_term) + 
		(viscosity_term * params.dynamic_viscosity);

	if(length(normal) > params.surface_tension_threshold) {
		sum += -params.surface_tension * color_field_laplacian * normal / length(normal);
	}
	
	return sum;
}