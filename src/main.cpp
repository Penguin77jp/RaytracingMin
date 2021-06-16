#include <iostream>
#include <string>
#include <filesystem>

#include "Renderer.h"
#include "SettingData.h"
#include "Camera.h"

int main(int argc, char* argv[]) {
	png::BiconvexLens lens(1,1,1);
	lens.GetRay(0,0);
	return 0;

	std::string jsonFile = "settingData.json";
	for (int i = 1; i < argc; ++i) {
		std::string command = argv[i];
		if (command == "-?") {
			std::cout << "Usage: RaytracingMin.exe [OPTION]..." << std::endl << std::endl
				<< "-sampleJson : サンプルシーンのjsonファイルを保存 ファイル名:settingData.json" << std::endl
				<< "-json : jsonファイルの指定 指定しない場合はsettingData.jsonで読み込む" << std::endl;
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
			std::cout << "サンプルシーンのjsonファイルを保存 ファイル名:settingData.json" << std::endl;
			png::LoadData::SaveSampleJson("settingData.json");
			return 0;
		}
	}

	if (!std::filesystem::exists(jsonFile)) {
		std::cout << "サンプルシーンのjsonファイルを保存 ファイル名:settingData.json" << std::endl;
		png::LoadData::SaveSampleJson(jsonFile);
	}
	std::cout << jsonFile << "を読み込みレンダリング開始" << std::endl;
	png::LoadData loadData(jsonFile);
	png::Renderer renderer(loadData.data);

	renderer.Render("result");
	
	return 0;
}