#ifndef _COLLISION_COLUMES_LOADER_H_
#define _COLLISION_COLUMES_LOADER_H_

#include <fstream>
#include "../util/pico_json/picojson.h"
#include "common/structures.h"

class collision_volumes_loader {
public:
	collision_volumes load_standard_json(std::ifstream& json_stream);

	collision_box read_box_json(picojson::object o);
	collision_sphere read_sphere_json(picojson::object o);
	collision_capsule read_capsule_json(picojson::object o);
};

#endif