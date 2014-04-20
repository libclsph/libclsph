libclsph
========

An OpenCL based GPU accelerated SPH fluid simulation library

Can I see it in action?
-----------------------

[Demo #1](https://www.youtube.com/watch?v=s7RinwmAOLU)    
[Demo #2](https://www.youtube.com/watch?v=uJznjnEIfrg)

Why?
-----------------------

Libclsph was created to explore the possibilty of using the power of OpenCL to speed up the simulation of SPH fluid mechanics.

Smoothed particle hydrodynamics is a fluid simulation technique that can be used to produce realistic simulations for animation,CGI or videogames.

How ?
-----------------------

Libclsph uses C++11 and OpenCL to run simulations. Results are exported using industry standards and can then be used in programs like [houdini](http://www.sidefx.com/). Simulation properties can be easily modified by the user.


Getting started
----------------

Be sure to visit the [wiki](https://github.com/libclsph/libclsph/wiki) for more information and detailed instructions on how to get started.

Libraries used
----------------
* [picojson](https://github.com/kazuho/picojson) is used to load simulation properties from files.
* [cereal](http://uscilab.github.io/cereal/) is used for all serialization needs.
* [partio](http://www.disneyanimation.com/technology/partio.html) is used to import/export particle data.
* [tinyobjloader](https://github.com/syoyo/tinyobjloader) is used to load geometry in obj format.
