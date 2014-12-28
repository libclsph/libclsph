#define EXIT_ON_CL_ERROR

#include <iostream>
#include <string>

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include "sph_simulation.h"

int main(int, char**) {
    sph_simulation simulation;

    try{
        simulation.load_settings(
            std::string("fluid_properties/water.json"),
            std::string("simulation_properties/default.json"));
    }
    catch (const std::exception& ex) {
        std::cerr << ex.what() << std::endl;
        std::exit(-1);
    }

    GLFWwindow* window;
    if (!glfwInit()) {
        return -1;
    }

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    window = glfwCreateWindow(640, 480, "SPH example", NULL, NULL);
    if (!window) {
        glfwTerminate();
        return -1;
    }

    glfwMakeContextCurrent(window);
    GLenum err = glewInit();
    if (GLEW_OK != err) {
        std::cerr << "glew error" << std::endl;
    }

    GLuint VertexArrayID;
    glGenVertexArrays(1, &VertexArrayID);
    glBindVertexArray(VertexArrayID);

    GLuint vertexbuffer;
    glGenBuffers(1, &vertexbuffer);

    simulation.pre_frame = [&] (particle* particles, const simulation_parameters& params, bool, profile_data&) {
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
        glBufferData(
            GL_ARRAY_BUFFER, 
            params.particles_count * sizeof(particle), 
            particles, 
            GL_STATIC_DRAW);

        glEnableVertexAttribArray(0);
        glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
        glVertexAttribPointer(
           0,
           3,
           GL_FLOAT,
           GL_FALSE,
           sizeof(particle),
           (void*)0);
        glDrawArrays(GL_POINTS, 0, params.particles_count);

        glDisableVertexAttribArray(0);

        glfwSwapBuffers(window);
        glfwPollEvents();
    };

    if(!simulation.current_scene.load("monkey.obj")) {
        std::cerr << "Unable to load scene: " << "monkey.obj" << std::endl;
        return -1;
    }

    simulation.simulate();

    glfwTerminate();
    return 0;
}
