#pragma once

#include "Ray.h"
#include "SceneObject.h"
#include "Camera.h"
#include "sceneLight.h"

#include <string>
#include <cmath>
#include <nlohmann/json.hpp>

namespace png {
	struct SettingData {
		int renderType, width, height, samples, superSamples, spectrumSamples;
		vec3 cameraOrigin, cameraTarget;
		double fov;
		WavePlane* wave;
		std::vector<SceneObject*> object;
		SceneLight* sceneLight;

		// external
		double cameraAperture;

		void AnimationUpdate(float time);
	};

	struct LoadData {
		LoadData(std::string);
		static void SaveSampleJson(std::string);

		static void to_json(nlohmann::json& json, const SettingData& data);
		static void from_json(const nlohmann::json& json, SettingData& data);

		SettingData data;
	};
}