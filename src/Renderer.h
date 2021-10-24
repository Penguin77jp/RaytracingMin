#pragma once

#include "SettingData.h"

#include <string>
#include <nlohmann/json.hpp>
#include <vector>

namespace png {
	class Renderer {
	public:
		Renderer(SettingData&);

		void Render(std::string);
	private:
		std::vector<double> image;
		SettingData& data;
	};
}