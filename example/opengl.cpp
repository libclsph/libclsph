#define EXIT_ON_CL_ERROR

#define SCREEN_RESOLUTION_X 1280
#define SCREEN_RESOLUTION_Y 720

#include <iostream>
#include <string>

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include "sph_simulation.h"
#include <stddef.h>

GLuint compile_shader(std::string path, GLuint shader_type){

    std::string shader_source = readKernelFile(path);
    int shader_source_length = shader_source.length();
    const char *shader_source_cstr = shader_source.c_str();

    GLuint shader = glCreateShader(shader_type);
    glShaderSource(shader, 1, &shader_source_cstr, &shader_source_length);
    glCompileShader(shader);

    GLint status;

    glGetShaderiv(shader, GL_COMPILE_STATUS, &status);
    if(status != GL_TRUE){

        char buffer[512];
        glGetShaderInfoLog(shader, 512, NULL, buffer);
        std::cout << std::string(buffer) << std::endl;

        throw std::runtime_error("Shader compilation failed!");
    }

    return shader;
}

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

    window = glfwCreateWindow(SCREEN_RESOLUTION_X, SCREEN_RESOLUTION_Y, "SPH example", NULL, NULL);
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

    //-----------------------------------------------------
    // Create shader program
    //-----------------------------------------------------
    GLuint vertexShader = compile_shader("shaders/vert.glsl",GL_VERTEX_SHADER);
    GLuint fragmentShader = compile_shader("shaders/frag.glsl",GL_FRAGMENT_SHADER);

    GLuint shaderProgram = glCreateProgram();
    glAttachShader(shaderProgram,vertexShader);
    glAttachShader(shaderProgram, fragmentShader);
    glLinkProgram(shaderProgram);
    glUseProgram(shaderProgram);

    //-----------------------------------------------------
    // Setup shader data
    //-----------------------------------------------------

    GLint base_density = glGetUniformLocation(shaderProgram, "base_density");
    glUniform1f(base_density,simulation.parameters.fluid_density);

    //-----------------------------------------------------
    // Create lambda for preframe processing
    //-----------------------------------------------------

    simulation.pre_frame = [&] (particle* particles, const simulation_parameters& params, bool) {
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
        glBufferData(
            GL_ARRAY_BUFFER,
            params.particles_count * sizeof(particle),
            particles,
            GL_DYNAMIC_DRAW);

        glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);

        // Specify the layout of the vertex data
        GLint posAttrib = glGetAttribLocation(shaderProgram, "position");
        glEnableVertexAttribArray(posAttrib);
        glVertexAttribPointer(
           posAttrib,
           3,
           GL_FLOAT,
           GL_FALSE,
           sizeof(particle),
           (void*)offsetof(particle,position));

        GLint accAttrib = glGetAttribLocation(shaderProgram, "density");
        glEnableVertexAttribArray(accAttrib);
        glVertexAttribPointer(accAttrib,
                              1,
                              GL_FLOAT,
                              GL_FALSE,
                              sizeof(particle),
                              (void*)(offsetof(particle,density)));


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
