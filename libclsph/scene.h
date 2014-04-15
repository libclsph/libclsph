#ifndef _SCENE_H_
#define _SCENE_H_

#include <vector>
#include <string>

class scene {
public:
	bool load(std::string filename);

	unsigned int face_count;
	std::vector<float> face_normals;
	std::vector<float> vertices;
	std::vector<unsigned int> indices;
};

#endif