libclsph
========

An OpenCL based GPU accelerated SPH fluid simulation library

![alt tag](http://i.imgur.com/iambaDc.png)

Why?
-----------------------

Libclsph was created to explore the possibilty of using the power of OpenCL to speed up the simulation of SPH fluid mechanics.

Smoothed particle hydrodynamics is a fluid simulation technique that can be used to produce realistic simulations for animation,CGI or videogames.

How ?
-----------------------

Libclsph uses C++11 and OpenCL to run simulations. Results are exported using industry standards and can then be used in programs like [houdini](http://www.sidefx.com/). Simulation properties can be easily modified by the user.


Getting started
-----------------------

An exemple is provided and can be run this way :

> ./example.out \<fluid_name\> \<simulation_properties_name\> \<scene_name\> \<frames_folder_prefix\>

The OpenCL runtime is required to use the library. If no GPU is found it will default to using the CPU.

More information
----------------

Be sure to visit the [wiki](https://github.com/libclsph/libclsph/wiki) for more information.

Libraries used
----------------
[picojson](https://github.com/kazuho/picojson) is used to load simulation properties from files.

[cerea](http://uscilab.github.io/cereal/) is used for all serialization needs.


