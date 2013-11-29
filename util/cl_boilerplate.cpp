#include "cl_boilerplate.h"

cl_int init_cl_single_device(cl::Context* context, std::vector<cl::Device>& device_array, std::string platform_hint, std::string device_hint, bool logging_enabled) {
	cl_int cl_error;
    std::vector<cl::Platform> all_platforms;
    cl_error = cl::Platform::get(&all_platforms);

    if(cl_error != CL_SUCCESS) return cl_error;

    if(all_platforms.size() == 0) {
        if(logging_enabled) std::cerr << "No platforms available" << std::endl;
        return CL_INVALID_PLATFORM;
    }

    cl::Platform default_platform = all_platforms[0];
    bool found = false;
    for(auto& p: all_platforms) {
    	if(p.getInfo<CL_PLATFORM_NAME>() == platform_hint) {
    		if(logging_enabled) std::cout << "Found platform hinted" << std::endl;

    		found = true;
    		default_platform = p;
    		break;
    	}
    }

    if(logging_enabled && !found) std::cout << "Did not find hinted platform, defaulting to first available" << std::endl; 
    if(logging_enabled) std::cout << "Using platform: " << default_platform.getInfo<CL_PLATFORM_NAME>() << std::endl;

    std::vector<cl::Device> all_devices;
    cl_error = default_platform.getDevices(CL_DEVICE_TYPE_ALL, &all_devices);
    if(cl_error != CL_SUCCESS) return cl_error;

    if(all_devices.size() == 0) {
        if(logging_enabled) std::cerr << "No device found" << std::endl;
        return CL_INVALID_DEVICE;
    }

    cl::Device default_device = all_devices[0];
    found = false;
    for(auto& d: all_devices) {
    	if(d.getInfo<CL_DEVICE_NAME>() == device_hint) {
    		if(logging_enabled) std::cout << "Found device hinted" << std::endl;

    		found = true;
    		default_device = d;
    		break;
    	}
    }

    if(logging_enabled && !found) std::cout << "Did not find hinted device, defaulting to first available" << std::endl; 
    if(logging_enabled) std::cout << "Using device: " << default_device.getInfo<CL_DEVICE_NAME>() << std::endl;

    device_array.push_back(default_device);

    *context = cl::Context(device_array, NULL, NULL, NULL, &cl_error);

    if(cl_error != CL_SUCCESS) return cl_error;

    if(logging_enabled) std::cout << "Context created" << std::endl;

    return CL_SUCCESS;
}

cl_int make_program(cl::Program* program, cl::Context& context, std::vector<cl::Device>& device_array, std::string program_source, bool logging_enabled, const char* options) {
    cl::Program::Sources source;

    std::pair<const char*, size_t> kernel_info;
    kernel_info.first = program_source.c_str();
    kernel_info.second = program_source.length();
    source.push_back(kernel_info);

    if(logging_enabled) std::cout << "source created" << std::endl;

    *program = cl::Program(context, source);
    if(logging_enabled) std::cout << "program object created but not built" << std::endl;

    if(program->build(device_array, options) != CL_SUCCESS) {
        if(logging_enabled) std::cerr << "Error building program" << std::endl;
        if(logging_enabled) std::cerr << program->getBuildInfo<CL_PROGRAM_BUILD_LOG>(device_array[0]) << std::endl;
        return CL_INVALID_PROGRAM;
    }

    if(logging_enabled) std::cout << "program built" << std::endl;
    if(logging_enabled) std::cout << program->getBuildInfo<CL_PROGRAM_BUILD_LOG>(device_array[0]) << std::endl;

    return CL_SUCCESS;
}

std::string readKernelFile(std::string kernelFileName){

    std::ifstream stream(kernelFileName);

    if(!stream)
    {
        std::cout<< kernelFileName << " file not found." << std::endl;
        exit(1);
    }

    std::stringstream buffer;
    buffer << stream.rdbuf();

    std::string programSource = buffer.str();

    if(programSource.length()==0)
    {
        std::cout<< kernelFileName << " is empty." << std::endl;
        exit(1);
    }

    return programSource;
}
