#include "Camera.h"
#include <iostream>

namespace png {
	Camera::Camera(SettingData& data)
		: org(vec3(json["camera"]["origin"][0], json["camera"]["origin"][1], json["camera"]["origin"][2]))
		, tar(vec3(json["camera"]["target"][0], json["camera"]["target"][1], json["camera"]["target"][2]))
		, up(vec3(0, 0, 1))
		, fov(json["camera"]["fov"])
	{}

	float Camera::GetFOV() const { return fov; }
}