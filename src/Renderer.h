#pragma once

#include "SettingData.h"

#include <string>
#include <nlohmann/json.hpp>
#include <vector>
#include <utility>

namespace png {
	class Renderer {
	public:
		Renderer(SettingData&);

		void Render(std::string);
	public:
		std::vector<double> image;
		SettingData& data;
		void AnimationUpdate(int frame, int allFrame, const clock_t benchmarkTime, Random& random);

		// wave
		std::vector<std::pair<double, std::pair<double, double>>> waveSetting;
	};
}