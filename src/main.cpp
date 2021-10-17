#include <iostream>
#include <string>
#include <filesystem>

#include "Renderer.h"
#include "SettingData.h"

int main(int argc, char* argv[]) {
	std::string jsonFile = "settingData.json";
	for (int i = 1; i < argc; ++i) {
		std::string command = argv[i];
		if (command == "-?") {
			std::cout << "Usage: RaytracingMin.exe [OPTION]..." << std::endl << std::endl
				<< "-sampleJson : save sample json file. file name is settingData.json" << std::endl
				<< "-json : specify json file. If you don't specify json file, automatically specify settingData.json." << std::endl;
			return 0;
		}
		else if (command == "-json") {
			if (argc <= i + 1) {
				jsonFile = "settingData.json";
			}
			else if (argv[i + 1][0] == '-') {
				jsonFile = "settingData.json";
			}
			else {
				jsonFile = argv[i + 1];
			}
		}
		else if (command == "-sampleJson") {
			std::cout << "export sample json file." << std::endl;
			png::LoadData::SaveSampleJson("settingData.json");
			return 0;
		}
    }

	if (!std::filesystem::exists(jsonFile)) {
		std::cout << "export sample json file." << std::endl;
		png::LoadData::SaveSampleJson(jsonFile);
	}
	std::cout << jsonFile << " is loaded" << std::endl;
    png::LoadData loadData(jsonFile);
	png::Renderer renderer(loadData.data);

	renderer.Render("result");
	
	return 0;
}