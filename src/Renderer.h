#pragma once
#include <string>
#include <nlohmann/json.hpp>
#include <vector>

#include "Camera.h"

namespace png {
	class Renderer {
	public:
		Renderer(nlohmann::json&);

		void Render(Camera&, std::string);
	private:
		std::vector<unsigned char> image;
	};
}