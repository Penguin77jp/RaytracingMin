#pragma once
#include <string>
#include <nlohmann/json.hpp>
#include "Ray.h"

namespace png {
	struct Material {
		vec3 color;
		float kd;
		vec3 emission;
	};
	struct Object {
		vec3 position;
		float size;
		Material material;
	};
	struct Camera {
		vec3 origin, target, upVec;
		float fov;
	};
	struct SettingData {
		int width, height, sample;
		Camera camera;
		std::vector<Object> object;
	};

	struct LoadData{
		LoadData(std::string);
		static void SaveSampleJson(std::string);

		static void to_json(nlohmann::json& json, const SettingData& data);
		static void from_json(const nlohmann::json& json, SettingData& data);

		SettingData data;
	};
}