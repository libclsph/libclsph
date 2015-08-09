#ifndef _HOUDINI_FILE_SAVER_H_
#define _HOUDINI_FILE_SAVER_H_

#include <string>

#include "common/structures.h"

class houdini_file_saver {
 public:
  houdini_file_saver(std::string frames_folder_prefix)
      : frames_folder_prefix(frames_folder_prefix), frame_count(0) {}

  int writeFrameToFile(particle* particles,
                       const simulation_parameters& parameters);

  std::string frames_folder_prefix;

 private:
  int frame_count;
};

#endif