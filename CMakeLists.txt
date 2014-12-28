cmake_minimum_required(VERSION 2.6)
project(libclsph)

# User provided build options
option(BUILD_EXAMPLES "Build examples" ON)

# Check if OpenCL is present on the system
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
find_package(OpenCL REQUIRED)
include_directories(${OPENCL_INCLUDE_DIRS})

# Set common compiler flags
set(COMPILE_FLAGS "--std=c++11 -g -Wfatal-errors -Wall -Wextra -Wno-unused-parameter -Wno-sign-compare -Wno-deprecated-declarations")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${COMPILE_FLAGS}")

# Copy headers
file(COPY libclsph/sph_simulation.h DESTINATION ${PROJECT_BINARY_DIR}/include/clsph)
file(COPY libclsph/scene.h DESTINATION ${PROJECT_BINARY_DIR}/include/clsph)
file(COPY libclsph/common DESTINATION ${PROJECT_BINARY_DIR}/include/clsph)
file(COPY libclsph/file_save_delegates/houdini_file_saver.h DESTINATION ${PROJECT_BINARY_DIR}/include/clsph/file_save_delegates)
file(COPY util DESTINATION ${PROJECT_BINARY_DIR}/include/clsph)

# Copy kernels and properties files
file(COPY libclsph/kernels DESTINATION ${PROJECT_BINARY_DIR}/examples)
file(COPY simulation_properties DESTINATION ${PROJECT_BINARY_DIR}/examples)
file(COPY scenes DESTINATION ${PROJECT_BINARY_DIR}/examples)
file(COPY fluid_properties DESTINATION ${PROJECT_BINARY_DIR}/examples)
file(COPY libclsph/common DESTINATION ${PROJECT_BINARY_DIR}/examples)

include_directories(${PROJECT_BINARY_DIR}/include/clsph)

# Build the main lib
add_library(clsph 
    libclsph/sph_simulation.cpp
    libclsph/file_save_delegates/houdini_file_saver.cpp
    util/cl_boilerplate.cpp
    util/houdini_geo/HoudiniFileDumpHelper.cpp
    libclsph/scene.cpp
    util/tinyobj/tiny_obj_loader.cc)

set_target_properties(clsph PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set_target_properties(clsph PROPERTIES ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)

target_link_libraries(clsph ${OPENCL_LIBRARIES})

if(BUILD_EXAMPLES)
    link_directories(${PROJECT_BINARY_DIR}/lib)

    # Houdini output example
    add_executable(particles example/particles.cpp)
    set_target_properties(particles PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/examples)
    target_link_libraries(particles clsph)

endif(BUILD_EXAMPLES)