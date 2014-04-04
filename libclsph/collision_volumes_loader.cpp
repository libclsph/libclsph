#include "collision_volumes_loader.h"

collision_volumes collision_volumes_loader::load_standard_json(std::ifstream& json_stream) {
	collision_volumes volumes;

	for(int i = 0; i < COLLISION_VOLUMES_COUNT; ++i) {
		volumes.spheres[i].active = 0;
		volumes.boxes[i].active = 0;
		volumes.capsules[i].active = 0;
	}

	picojson::value v;
	json_stream >> v;

	picojson::object volumes_object = v.get<picojson::object>();

	volumes.scene_bounding_box = this->read_box_json(volumes_object["bounding_box"].get<picojson::object>());

	int index = 0;
	const picojson::array& spheres = volumes_object["spheres"].get<picojson::array>();
    for (picojson::array::const_iterator i = spheres.begin(); i != spheres.end(); ++i) {
    	volumes.spheres[index] = this->read_sphere_json((*i).get<picojson::object>());
        ++index;
    }

    index = 0;
     picojson::array& boxes = volumes_object["boxes"].get<picojson::array>();
    for (picojson::array::const_iterator i = boxes.begin(); i != boxes.end(); ++i) {
        volumes.boxes[index] = this->read_box_json((*i).get<picojson::object>());
        ++index;
    }

    index = 0;
    const picojson::array& capsules = volumes_object["capsules"].get<picojson::array>();
    for (picojson::array::const_iterator i = capsules.begin(); i != capsules.end(); ++i) {
        volumes.capsules[index] = this->read_capsule_json((*i).get<picojson::object>());
        ++index;
    }

	return volumes;
}

collision_box collision_volumes_loader::read_box_json(picojson::object o) {
	collision_box b;

    b.center.s[0] = (float)o["center"]
        .get<picojson::object>()["x"]
        .get<double>();
    b.center.s[1] = (float)o["center"]
        .get<picojson::object>()["y"]
        .get<double>();
    b.center.s[2] = (float)o["center"]
        .get<picojson::object>()["z"]
        .get<double>();

    b.axis_extends.s[0] = (float)o["axis_extends"]
        .get<picojson::object>()["x"]
        .get<double>();
    b.axis_extends.s[1] = (float)o["axis_extends"]
        .get<picojson::object>()["y"]
        .get<double>();
    b.axis_extends.s[2] = (float)o["axis_extends"]
        .get<picojson::object>()["z"]
        .get<double>();

    b.container_or_obstacle = o["container"].get<bool>() ? 1 : -1;
    b.active = 1;

    return b;
}

collision_sphere collision_volumes_loader::read_sphere_json(picojson::object o) {
	collision_sphere s;

    s.center.s[0] = (float)o["center"]
        .get<picojson::object>()["x"]
        .get<double>();
    s.center.s[1] = (float)o["center"]
        .get<picojson::object>()["y"]
        .get<double>();
    s.center.s[2] = (float)o["center"]
        .get<picojson::object>()["z"]
        .get<double>();

    s.radius = (float)o["radius"]
        .get<double>();

    s.container_or_obstacle = o["container"].get<bool>() ? 1 : -1;
    s.active = 1;

    return s;
}

collision_capsule collision_volumes_loader::read_capsule_json(picojson::object o) {
	collision_capsule c;

    c.p0.s[0] = (float)o["p0"]
        .get<picojson::object>()["x"]
        .get<double>();
    c.p0.s[1] = (float)o["p0"]
        .get<picojson::object>()["y"]
        .get<double>();
    c.p0.s[2] = (float)o["p0"]
        .get<picojson::object>()["z"]
        .get<double>();

    c.p1.s[0] = (float)o["p1"]
        .get<picojson::object>()["x"]
        .get<double>();
    c.p1.s[1] = (float)o["p1"]
        .get<picojson::object>()["y"]
        .get<double>();
    c.p1.s[2] = (float)o["p1"]
        .get<picojson::object>()["z"]
        .get<double>();

    c.radius = (float)o["radius"]
        .get<double>();

    c.container_or_obstacle = o["container"].get<bool>() ? 1 : -1;
    c.active = 1;

	return c;
}