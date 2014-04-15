#include "scene.h"

#include "../util/tinyobj/tiny_obj_loader.h"

#include <iostream>
#include <cassert>
#include <cmath>

bool scene::load(std::string filename) {
	std::vector<tinyobj::shape_t> shapes;

	std::string err = tinyobj::LoadObj(shapes, (std::string("scenes/") + filename).c_str(), "scenes/");

	if(!err.empty()) {
		std::cerr << err << std::endl;
		return false;
	}

	std::cout << 
		"Scene Loading - number of shapes in file [" << filename << "]: " << shapes.size() << std::endl;

	for(size_t i = 0; i < shapes.size(); ++i) {

		//TODO: handle multiple shapes properly
		indices = shapes[i].mesh.indices;
		if(indices.size() % 3 != 0) {
			std::cerr << "Meshes must be made of triangles only" << std::endl;
			return false;
		}

		face_count = indices.size() / 3;

		vertices = shapes[i].mesh.positions;

		//Pre-compute face normals
		for(size_t j = 0; j < face_count; ++j) {

			size_t off = 3 * j;
			float nx, ny, nz;
			float ux, uy, uz, vx, vy, vz;

			ux = vertices[3 * indices[off + 1] + 0] - vertices[3 * indices[off + 0] + 0];
			uy = vertices[3 * indices[off + 1] + 1] - vertices[3 * indices[off + 0] + 1];
			uz = vertices[3 * indices[off + 1] + 2] - vertices[3 * indices[off + 0] + 2];

			vx = vertices[3 * indices[off + 2] + 0] - vertices[3 * indices[off + 0] + 0];
			vy = vertices[3 * indices[off + 2] + 1] - vertices[3 * indices[off + 0] + 1];
			vz = vertices[3 * indices[off + 2] + 2] - vertices[3 * indices[off + 0] + 2];

			nx = uy * vz - uz * vy;
			ny = uz * vx - ux * vz;
			nz = ux * vy - uy * vx;

			float length = sqrt(nx * nx + ny * ny + nz * nz);

			face_normals.push_back(nx / length);
			face_normals.push_back(ny / length);
			face_normals.push_back(nz / length);
		}
	}
	return true;
}