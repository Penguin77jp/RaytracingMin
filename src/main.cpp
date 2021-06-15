#include <iostream>
#include <string>

#include "Renderer.h"
#include "SettingData.h"

int main(int argc, char* argv[]) {
	std::string jsonFile = "settingData.json";
	for (int i = 1; i < argc; ++i) {
		std::string command = argv[i];
		std::cout << command << std::endl;
		if (command == "-?") {
			std::cout << "Usage: RaytracingMin.exe [OPTION]..." << std::endl
				<< "-sampleJson : �T���v���V�[����json�t�@�C����ۑ� �t�@�C����:settingData.json" << std::endl
				<< "-json : json�t�@�C���̎w�� �w�肵�Ȃ��ꍇ��settingData.json�œǂݍ���" << std::endl;
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
			std::cout << "�T���v���V�[����json�t�@�C����ۑ� �t�@�C����:settingData.json" << std::endl;
			png::LoadData::SaveSampleJson("settingData.json");
			return 0;
		}
	}

	std::cout << jsonFile << "��ǂݍ��݃����_�����O�J�n" << std::endl;
	png::LoadData loadData(jsonFile);
	png::Renderer renderer(loadData.data);

	renderer.Render("result");
	
	return 0;
}