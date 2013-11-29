#ifndef _CL_BOILERPLATE_H_
#define _CL_BOILERPLATE_H_

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>

#ifdef __APPLE__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-pedantic"
#pragma GCC diagnostic ignored "-Wall"
#pragma GCC diagnostic ignored "-Wextra"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#include "../util/cl.hpp"
#pragma GCC diagnostic pop
#else
#include <CL/cl.hpp>
#endif

#ifdef EXIT_ON_CL_ERROR
#define check_cl_error exit_on_error
#else
#define check_cl_error check_error
#endif

#define check_error(error) error == CL_SUCCESS ? true : \
    ([=](cl_int err) \
        { std::cerr << "An OpenCL error occured (" << __FILE__ << ":" << __LINE__ << ")-> " << err << std::endl; return false; })(error);

#define exit_on_error(error) error == CL_SUCCESS ? true : \
    ([=](cl_int err) \
        { std::cerr << "An OpenCL error occured (" << __FILE__ << ":" << __LINE__ << ")-> " << err << std::endl; exit(-1); return false; })(error);

/*
	Initializes OpenCL and creates a context with one device attached.  Also populates device_array with the device used
	Attempts to select the platform and device that have platform_hint and device_hint for names.
	setting logging_enabled to true will cause debug and error info to be logged to std::cout and std::cerr
*/
cl_int init_cl_single_device(cl::Context* context, std::vector<cl::Device>& device_array, std::string platform_hint = "", std::string device_hint = "", bool logging_enabled = false);
cl_int make_program(cl::Program* program, cl::Context& context, std::vector<cl::Device>& device_array, std::string program_source, bool logging_enabled = false, const char* options = "");
std::string readKernelFile(std::string kernelFileName);

#endif