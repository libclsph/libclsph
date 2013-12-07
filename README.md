libclsph
========

An OpenCL based GPU accelerated SPH fluid simulation library



How to use the library:
-----------------------

An exemple is provided and can be run this way :

> ./example.out \<fluid_name\> \<simulation_properties_name\> \<scene_name\> \<frames_folder_prefix\>

The OpenCL runtime is required to use the library. If no GPU is found it will default to using the CPU.
