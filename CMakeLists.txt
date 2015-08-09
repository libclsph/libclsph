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

# Copy kernels, shader and properties files
file(COPY libclsph/kernels DESTINATION ${PROJECT_BINARY_DIR}/examples)
file(COPY shaders DESTINATION ${PROJECT_BINARY_DIR}/examples)
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
    add_executable(clsphparticles example/particles.cpp)
    set_target_properties(clsphparticles PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/examples)
    target_link_libraries(clsphparticles clsph)

    find_package(OpenGL)
    find_package(GLEW)
    find_package(GLFW)

    # OpenGL rendered example requires OpenGL, GLEW and GLFW
    if(OPENGL_FOUND AND GLEW_FOUND AND GLFW_FOUND)
        include_directories(${OPENGL_INCLUDE_DIRS})
        include_directories(${GLEW_INCLUDE_DIRS})
        include_directories(${GLFW_INCLUDE_DIRS})

        add_executable(clsphopengl example/opengl.cpp)
        set_target_properties(clsphopengl PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/examples)

        target_link_libraries(clsphopengl clsph)
        target_link_libraries(clsph ${OPENGL_LIBRARIES})
        target_link_libraries(clsph ${GLEW_LIBRARIES})
        target_link_libraries(clsph ${GLFW_LIBRARIES} -lX11 -lXi -lXrandr -lXxf86vm -lXinerama -lXcursor -lrt -lm -pthread )
    endif(OPENGL_FOUND AND GLEW_FOUND AND GLFW_FOUND)
endif(BUILD_EXAMPLES)
