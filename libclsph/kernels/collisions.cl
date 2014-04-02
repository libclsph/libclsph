#include "common.cl"

typedef struct {
	float3 position, next_velocity;
} collision_response;

typedef struct {
	float3 collision_point, surface_normal;
	float penetration_depth;
} collision;

float sphere(float3 point, collision_sphere s);
float capsule(float3 point, collision_capsule c);
float box(float3 point, collision_box b);

collision interpret_sphere_collision(float3 position, float func, collision_sphere s);
collision interpret_capsule_collision(float3 position, float func, collision_capsule c);
collision interpret_box_collision(float3 position, float func, collision_box b);

int respond(collision_response* response, collision c, simulation_parameters params);
collision_response handle_collisions(float3 position, float3 next, simulation_parameters params, collision_volumes volumes);

float sphere(float3 point, collision_sphere s) {
	return pown(length(point - s.center), 2) - pown(s.radius, 2);
}

float capsule(float3 point, collision_capsule c) {
	float3 q = c.p0 + fmin(1.f, fmax(0.f, -dot(c.p0 - point, c.p1 - c.p0) / pown(length(c.p1 - c.p0), 2))) * (c.p1 - c.p0);
	return length(q - point) - c.radius;
}

inline float max_c(float3 v) {
	return fmax(v.x, fmax(v.y, v.z));
}

float box(float3 point, collision_box b) {
	//TODO: allow orientations and project point into obb space (see kelager)
	return max_c(fabs(point - b.center) - b.axis_extends);
}

collision interpret_sphere_collision(float3 position, float func, collision_sphere s) {
	collision c;

	c.collision_point = 
		s.center + 
		s.radius *
		(position - s.center) /
		length(position - s.center);

	c.penetration_depth = length(s.center - position);
	c.surface_normal = 
		sign(func) * 
		(s.center - position) / 
		length(s.center - position);

	return c;
}

collision interpret_capsule_collision(float3 position, float func, collision_capsule c) {
	collision col;

	//TODO: compute only once
	float3 q = c.p0 + fmin(1.f, fmax(0.f, -dot(c.p0 - position, c.p1 - c.p0) / pown(length(c.p1 - c.p0), 2))) * (c.p1 - c.p0);

	col.collision_point = 
		q + c.radius * (position - q) / length(position - q);

	col.penetration_depth = fabs(capsule(position, c));
	col.surface_normal = sign(func) * (q - position) / length(q - position);

	return col;
}

collision interpret_box_collision(float3 position, float func, collision_box b) {
	collision c;

	float3 local_collision_point = fmin(b.axis_extends, fmax(-b.axis_extends, position - b.center));
	c.collision_point = b.center + local_collision_point;
		
	c.penetration_depth = length(c.collision_point - position);

	float3 s = sign(local_collision_point - (position - b.center));

	c.surface_normal = normalize(s);

	return c;
}

int respond(collision_response* response, collision c, simulation_parameters params) {
	response->position = c.collision_point;

	response->next_velocity -= 
		(1.f + 
			params.restitution * 
			c.penetration_depth / 
			(params.time_delta * params.simulation_scale * length(response->next_velocity))) * 
		dot(response->next_velocity, c.surface_normal) * c.surface_normal;

	//response->next_velocity *= 0.75f;

	return 1;
}

collision_response handle_collisions(float3 position, float3 next, simulation_parameters params, collision_volumes volumes) {
	collision_response response = {
		position, next,
	};

	for(int i = 0; i < COLLISION_VOLUMES_COUNT; ++i) {
		if(volumes.spheres[i].active) {
			float s = sphere(response.position, volumes.spheres[i]);
			if((s * volumes.spheres[i].container_or_obstacle) > 0.f) {
				collision col = interpret_sphere_collision(response.position, s, volumes.spheres[i]);
				respond(&response, col, params);
			}
		}

		if(volumes.capsules[i].active) {
			float c = capsule(response.position, volumes.capsules[i]);
			if((c * volumes.capsules[i].container_or_obstacle) > 0.f) {
				collision col = interpret_capsule_collision(response.position, c, volumes.capsules[i]);
				respond(&response, col, params);
			}
		}

		if(volumes.boxes[i].active) {
			float b = box(response.position, volumes.boxes[i]);
			if(b * volumes.boxes[i].container_or_obstacle > 0.f) {
				collision col = interpret_box_collision(response.position, b, volumes.boxes[i]);
				respond(&response, col, params);
			}
		}
	}

	float b = box(response.position, volumes.scene_bounding_box);
	if(b * volumes.scene_bounding_box.container_or_obstacle > 0.f) {
		collision col = interpret_box_collision(response.position, b, volumes.scene_bounding_box);
		respond(&response, col, params);
	}

	return response;
}