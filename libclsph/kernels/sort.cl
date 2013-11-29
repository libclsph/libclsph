constant const unsigned int mask = 0xFF;

inline unsigned int get_count_offset(int index, unsigned int mask, int pass_number, int radix_width) {
	return (index & (mask << (pass_number * radix_width))) >> (pass_number * radix_width);
}

inline uint2 get_start_and_end(size_t particle_count, int thread_count, int work_item_id) {

	size_t particles_per_thread = particle_count / thread_count;
	size_t start_index = particles_per_thread * work_item_id;
	size_t end_index = start_index + particles_per_thread - 1;

	if(work_item_id == thread_count - 1) {
		end_index = particle_count - 1;
	}

	return (uint2)(start_index, end_index);
}

/* size of counts = sizeof(size_t) * bucket_count * thread_count */
void kernel sort_count(
	global const particle* particles, 
	global volatile unsigned int* counts, 
	simulation_parameters params,
	int thread_count,
	int pass_number,
	int radix_width) {

	size_t work_item_id = get_global_id(0);
	uint2 indices = get_start_and_end(params.particles_count, thread_count, work_item_id);

	for(size_t i = indices.x; i <= indices.y; ++i) {
		unsigned int bucket = get_count_offset(particles[i].grid_index, mask, pass_number, radix_width);

		size_t counts_index = bucket * thread_count + work_item_id;

		++(counts[counts_index]);
	}
}

void kernel sort(
	global const particle* in_particles, 
	global particle* out_particles, 
	global size_t* start_indices, 
	simulation_parameters params,
	int thread_count,
	int pass_number,
	int radix_width) {
	
	size_t work_item_id = get_global_id(0);
	uint2 indices = get_start_and_end(params.particles_count, thread_count, work_item_id);

	for(size_t i = indices.x; i <= indices.y; ++i) {
		unsigned int bucket = get_count_offset(in_particles[i].grid_index, mask, pass_number, radix_width);

		size_t insertion_index = start_indices[bucket * thread_count + work_item_id];
		++(start_indices[bucket * thread_count + work_item_id]);

		out_particles[insertion_index] = in_particles[i];
	}
}