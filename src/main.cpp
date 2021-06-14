#include <nlohmann/json.hpp>
#include <fstream>

#include "Camera.h"
#include "Renderer.h"

int main() {
	std::ifstream jsonStream("settingData.json");
	nlohmann::json jsonSetting;
	jsonStream >> jsonSetting;
	png::Camera cam(jsonSetting);
	png::Renderer renderer(jsonSetting);
	renderer.Render(cam,"result.png");
	
	return 0;
}