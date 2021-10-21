#pragma once
#include <string>
#include <cmath>
#include <nlohmann/json.hpp>
#include "Ray.h"
#include "SceneObject.h"

namespace png {
	struct Camera {
		vec3 origin, target, upVec;
		float fov;
	};
	struct SettingData {
		int width, height, samples, superSamples, spectrumSamples;
		Camera camera;
		std::vector<SceneObject*> object;
	};

	struct LoadData {
		LoadData(std::string);
		static void SaveSampleJson(std::string);

		static void to_json(nlohmann::json& json, const SettingData& data);
		static void from_json(const nlohmann::json& json, SettingData& data);

		SettingData data;
	};
}