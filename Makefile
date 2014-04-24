#Target platform string.  Either osx or linux
PLATFORM=linux

#Path to OpenCL headers' containing folder
OPENCL_INCLUDE_DIR=/opt/AMDAPP/include 
#OPENCL_INCLUDE_DIR=/usr/local/NVIDIA_OpenCL_SDK/src/NVIDIA\ GPU\ Computing\ SDK/OpenCL/common/inc/

#Path to OpenCL library's containing folder
#OPENCL_LIB_DIR=/opt/AMDAPP/lib/x86_64
OPENCL_LIB_DIR=/usr/local/lib 

#Desired output directory of the build
OUT_DIR = bin/

#Path to partio, leave empty to build without
LIB_PARTIO_A_PATH=

#Uncomment this line if a local cl.hpp is required
#LOCAL_CL_HPP=-DUSE_LOCAL_CL_HPP

#No editing should be required beyond this line

#If a path for partio is provided then define USE_PARTIO
ifneq ($(LIB_PARTIO_A_PATH),)
  DEFINES=-DUSE_PARTIO
endif

INCLUDE_PATH = -I$(OPENCL_INCLUDE_DIR)
LIB_PATH = -L$(OPENCL_LIB_DIR)
ifeq ($(PLATFORM),osx)
	LIBS=-framework OpenCL -lz
else
	LIBS=-lOpenCL -lz
endif

FLAGS = --std=c++11 -g -Wfatal-errors -Wall -Wextra -Wno-sign-compare -Wno-deprecated-declarations $(LOCAL_CL_HPP)

LIBRARY_FILES = 											\
	libclsph/sph_simulation.cpp 							\
	libclsph/file_save_delegates/houdini_file_saver.cpp 	\
	util/cl_boilerplate.cpp 								\
	util/houdini_geo/HoudiniFileDumpHelper.cpp				\
	libclsph/scene.cpp 										\
	util/tinyobj/tiny_obj_loader.cc


_OBJECT_FILES =					\
	sph_simulation.o			\
	houdini_file_saver.o 		\
	cl_boilerplate.o			\
	HoudiniFileDumpHelper.o		\
	scene.o 					\
	tiny_obj_loader.o 			

PROP_FOLDERS = fluid_properties simulation_properties scenes
KERNEL_FOLDERS = libclsph/kernels libclsph/common

OBJECT_FILES = $(patsubst %, $(OUT_DIR)%, $(_OBJECT_FILES))

all: example/particles.cpp lib
	g++ $(INCLUDE_PATH) -I./libclsph $(LIB_PATH) $(FLAGS) example/particles.cpp -L./bin -lclsph $(LIBS) $(DEFINES) $(LIB_PARTIO_A_PATH) -o $(OUT_DIR)example.out
	cp -r $(PROP_FOLDERS) $(OUT_DIR)

lib: $(LIBRARY_FILES)
	g++ -c $(INCLUDE_PATH) $(LIB_PATH) $(FLAGS)  $(LIBRARY_FILES) $(DEFINES) $(LIB_PARTIO_A_PATH) $(LIBS)
	mkdir -p $(OUT_DIR)
	mv $(_OBJECT_FILES) $(OUT_DIR)
	ar rcs $(OUT_DIR)libclsph.a $(OBJECT_FILES)
	cp -r $(KERNEL_FOLDERS) $(OUT_DIR)

clean: 
	rm $(OBJECT_FILES)

clobber: clean
	rm $(OUT_DIR)libclsph.a $(OUT_DIR)example.out
