#pragma once
#include <string>
#include <cmath>
#include <nlohmann/json.hpp>
#include "Ray.h"

namespace png {
	struct Material {
		vec3 color;
		vec3 emission;
		vec3 colorKD() const {
			return color / kd();
		}
		float kd() const {
			return std::max(std::max(color.x, color.y), color.z);
		}
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
		int width, height, samples, superSamples;
		Camera camera;
		std::vector<Object> object;
	};

	struct LoadData {
		LoadData(std::string);
		static void SaveSampleJson(std::string);

		static void to_json(nlohmann::json& json, const SettingData& data);
		static void from_json(const nlohmann::json& json, SettingData& data);

		SettingData data;
	};
}