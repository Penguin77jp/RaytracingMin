#pragma once

#include "Ray.h"
#include "SceneObject.h"
#include "Camera.h"

#include <string>
#include <cmath>
#include <nlohmann/json.hpp>
#include <tinyxml2.h>

namespace png {
	struct SettingData {
		int renderType, width, height, samples, superSamples, spectrumSamples;
		vec3 cameraOrigin, cameraTarget;
		double fov;
		std::vector<SceneObject*> object;
	};

	struct LoadData {
		LoadData(std::string);
		//static void SaveSampleJson(std::string);

		void to_xml(const SettingData& in_data, tinyxml2::XMLDocument& out_doc);
		void from_xml(tinyxml2::XMLDocument& in_doc, SettingData& out_data);
		void Print();

		SettingData data;
	};
}