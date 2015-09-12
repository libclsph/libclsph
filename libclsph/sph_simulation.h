#ifndef _SPH_SIMULATION_H_
#define _SPH_SIMULATION_H_

#include <functional>
#include "common/structures.h"
#include "scene.h"

class sph_simulation {
 public:
  sph_simulation() : parameters(), write_intermediate_frames(false) {}

  void simulate(int frame_count = 0);

  simulation_parameters parameters;
  precomputed_kernel_values precomputed_terms;

  std::function<void(particle *, const simulation_parameters &, bool)>
      pre_frame;
  std::function<void(particle *, const simulation_parameters &, bool)>
      post_frame;

  void load_settings(std::string fluid_file_name,
                     std::string parameters_file_name);

  bool write_intermediate_frames;
  bool serialize;
  float initial_volume;
  scene current_scene;

 private:
  void init_particles(particle *buffer, const simulation_parameters &);
  void sort_particles(particle *, cl::Buffer &, cl::Buffer &, unsigned int *);
  void simulate_single_frame(particle *, particle *);

  cl::Context context_;
  cl::CommandQueue queue_;

  cl::Kernel kernel_density_pressure_;
  cl::Kernel kernel_advection_collision_;
  cl::Kernel kernel_forces_;
  cl::Kernel kernel_locate_in_grid_;
  cl::Kernel kernel_sort_count_;
  cl::Kernel kernel_sort_;

  cl::Buffer front_buffer_;
  cl::Buffer back_buffer_;

  static const int kSortThreadCount = 128;
  static const int kBucketCount = 256;
  static const int kRadixWidth = 8;

  std::array<unsigned int, kSortThreadCount * kBucketCount> sort_count_array_;
  cl::Buffer sort_count_buffer_;
};

#endif
