#pragma once
#include <string>
#include <nlohmann/json.hpp>
#include <vector>

#include "SettingData.h"

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