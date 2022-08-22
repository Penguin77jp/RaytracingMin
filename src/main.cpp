#include <iostream>
#include <string>
#include <filesystem>
#include <omp.h>

#include "Renderer.h"
#include "SettingData.h"

//debug
#include <iomanip>
#include <chrono>
#include "Color.h"
#include "Spectrum.h"
#include "WaveSolver.h"
#include <cmath>
#include "stb_image_write.h"
#include <format>

int main(int argc, char* argv[]) {
	/*
	WaveSolver wave;
	auto resultImage = std::vector<unsigned char>(wave.widthN() * wave.heightN() * 3);
	for (int f = 0; f < 10; ++f) {
		wave.step(wave.time2step(2.0));
		for (int y = 1; y < wave.heightN()-1; ++y) {
			for (int x = 1; x < wave.widthN()-1; ++x) {
				auto imageIndex = y * wave.widthN() * 3 + x * 3;
				auto imageVal = wave.normal(x, y) * 0.5 + png::vec3(0.5, 0.5, 0.5);
				auto imageValCharR = (unsigned char)(255.0 * std::min(std::max(imageVal.x, 0.0), 1.0));
				auto imageValCharG = (unsigned char)(255.0 * std::min(std::max(imageVal.y, 0.0), 1.0));
				auto imageValCharB = (unsigned char)(255.0 * std::min(std::max(imageVal.z, 0.0), 1.0));
				resultImage[imageIndex] = imageValCharR;
				resultImage[imageIndex + 1] = imageValCharG;
				resultImage[imageIndex + 2] = imageValCharB;
			}
			//stbi_write_jpg(std::format("wave{}.jpg", f).c_str(), wave.widthN(), wave.heightN(), 3, resultImage.data(), 100);
			stbi_write_bmp(std::format("wave{}.bmp", f).c_str(), wave.widthN(), wave.heightN(), 3, resultImage.data());
		}
	}
	*/

#ifdef _DEBUG
	std::cout << "==============DEBUG MODE==============" << std::endl;
#endif
	// max threads
	std::cout << "max threads : " << omp_get_max_threads() << std::endl;

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
	std::cout << jsonFile << " loaded" << std::endl;
	png::LoadData loadData(jsonFile);
	png::Renderer renderer(loadData.data);
	for (int f = 0; f < 50; ++f) {
		auto& obj = renderer.data.object;
		for (int i = 0; i < obj.size(); ++i) {
			obj[i]->AnimationUpdate(1);
		}
		std::cout << std::format("{} / {}", f, 50) << std::endl;
		renderer.Render(std::format("wave{}", f));
		renderer.image = std::vector<double>(renderer.image.size(), 0);
	}

	return 0;
}