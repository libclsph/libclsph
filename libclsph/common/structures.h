#ifndef _STRUCTURES_H_
#define _STRUCTURES_H_

#ifdef KERNEL_INCLUDE
#define cl_float float
#define cl_uint unsigned int
#define cl_int int
#define cl_float3 float3
#else
#include "../../util/cl_boilerplate.h"
#endif

#define COLLISION_VOLUMES_COUNT 3

typedef struct {
    cl_uint particles_count;
    cl_float max_velocity;
    cl_float fluid_density;
    cl_float total_mass;
    cl_float particle_mass;
    cl_float dynamic_viscosity;
    cl_float simulation_time;
    cl_float target_fps;
    cl_float h;
    cl_float simulation_scale;
    cl_float time_delta;
    cl_float surface_tension_threshold;
    cl_float surface_tension;
    cl_float restitution;
    cl_float K;

    cl_float3 constant_acceleration;

    cl_int grid_size_x;
    cl_int grid_size_y;
    cl_int grid_size_z;
    cl_uint grid_cell_count;
    cl_float3 min_point, max_point;
} simulation_parameters;

typedef struct {
    cl_float3 position, velocity, intermediate_velocity, constant_acceleration;
    cl_float density, pressure;
    cl_int grid_index;
} particle;

typedef struct {
    cl_float3 center;
    cl_float radius;
    cl_int container_or_obstacle; //-1 for obstacle, 1 for container
    cl_int active;
} collision_sphere;

typedef struct {
    cl_float3 p0;
    cl_float3 p1;
    cl_float radius;
    cl_int container_or_obstacle;
    cl_int active;
} collision_capsule;

typedef struct {
    cl_float3 center;
    cl_float3 axis_extends;
    cl_int container_or_obstacle;
    cl_int active;
} collision_box;

typedef struct {
    collision_box scene_bounding_box;
    collision_box boxes[COLLISION_VOLUMES_COUNT];
    collision_sphere spheres[COLLISION_VOLUMES_COUNT];
    collision_capsule capsules[COLLISION_VOLUMES_COUNT];
} collision_volumes;

#endif