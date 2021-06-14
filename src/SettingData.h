#pragma once
#include <string>
#include "Ray.h"

namespace png {
	struct SettingData{
		SettingData(std::string);

		int width, height, sample;
		vec3 org, tar, up;
		float fov;
	};
}