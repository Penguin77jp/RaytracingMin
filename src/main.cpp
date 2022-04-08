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

#include <stdio.h>
#include <vector>
#include <cmath> //for sin, cos

#include <iostream>
#include <vector>
#include <optional>
#include <fstream>
#include "Ray.h"
#include "2D.h"

struct Lens {
	double r, thickness, ior, apertureDiameter, z;
};

void Intersect(png::vec3 sphereCenter, Lens& lensData, png::Ray ray, std::optional<png::vec3>& hitpoint, png::vec3& normal) {
	if (std::abs(lensData.r) > 0.0) {
		auto p_o = sphereCenter - ray.org;
		auto b = Dot(p_o, ray.dir);
		double d4 = b * b - Dot(p_o, p_o) + pow2(lensData.r);

		if (d4 < 0.0) {
			hitpoint = std::nullopt;
			return;
		}

		double d4sqrt = std::sqrt(d4);
		double t1 = b - d4sqrt;
		double t2 = b + d4sqrt;

		double t = 0.0;
		if (t1 > png::EPS) {
			t = t1;
		}
		else if (t2 > png::EPS) {
			t = t2;
		}
		else {
			hitpoint = std::nullopt;
			return;
		}

		hitpoint = ray.org + ray.dir * t;
		if (pow2(hitpoint.value().x) + pow2(hitpoint.value().y) > pow2(lensData.apertureDiameter)) {
			hitpoint = std::nullopt;
			return;
		}

		normal = Normalize(hitpoint.value() - sphereCenter);
		if (Dot(ray.dir, normal) > 0) {
			normal = -normal;
		}
	}
	else {
		double t = (lensData.z - ray.org.z) / ray.dir.z;
		auto point = ray.org + t * ray.dir;
		if (pow2(hitpoint.value().x) + pow2(hitpoint.value().y) > pow2(lensData.apertureDiameter)) {
			hitpoint = std::nullopt;
			return;
		}

		normal = png::vec3(0, 0, -1);
		if (Dot(ray.dir, normal) > 0) {
			normal = -normal;
		}
	}
}

std::optional<png::vec3> Refract(png::vec3 incomingDir, png::vec3 lensNormal, double n1, double n2) {
	auto cos2 = 1.0 - pow2(n1 / n2) * (1.0 - pow2(Dot(-incomingDir, lensNormal)));
	if (cos2 <= 0.0) {
		// full refract
		return std::nullopt;
	}
	return n1 / n2 * incomingDir - (n1 / n2 * Dot(incomingDir, lensNormal) + std::sqrt(cos2)) * lensNormal;
}

std::vector<png::vec3> Pathtracing(std::vector<Lens>& lensData, png::Ray ray) {
	std::vector<png::vec3> points;
	double n1 = 1.0;
	points.push_back(ray.org);
	for (int i = 0; i < lensData.size(); ++i) {
		auto& targetLens = lensData[i];

		std::optional<png::vec3> hitpoint;
		png::vec3 lensNormal;
		Intersect(png::vec3(0, 0, targetLens.z + targetLens.r), targetLens, ray, hitpoint, lensNormal);
		if (!hitpoint.has_value()) {
			return std::vector<png::vec3>();
		}
		points.push_back(hitpoint.value());

		auto refractedDir = Refract(ray.dir, lensNormal, n1, targetLens.ior);

		if (!refractedDir.has_value()) {
			return std::vector<png::vec3>();
		}
		ray.dir = refractedDir.value();
		ray.org = hitpoint.value();
		n1 = targetLens.ior;
	}
	points.push_back(ray.org + ray.dir);

	return points;
}

void Output(std::ofstream& file, std::vector<png::vec3>& points) {
	if (points.size() == 0) {
		return;
	}
	for (int i = 0; i < points.size() - 1; ++i) {
		file << points[i].x << " " << points[i].y << " " << points[i].z << std::endl;
		file << points[i + 1].x << " " << points[i + 1].y << " " << points[i + 1].z << std::endl;
	}
}

/*
int main() {
	auto data = png::SettingData();
	data.width = 5;
	data.height = 5;
	png::LensSystem lensSystem("camera/convex.json", data);

	png::Ray ray;
	png::Ray incomingImageSensor;
	std::optional<png::Ray> generatedRay;
	png::Random rand;
	if (0 == 0) {
		for (double z = 0; z >= -100; z -= 10.0) {
			std::ofstream file("path" + std::to_string(int(z)) + ".txt");
			for (int i = 0; i < 100; ++i) {
				lensSystem.Test(2, 2, z, -1, incomingImageSensor, generatedRay, rand, file);
			}
		}
	}
	else {
		for (int y = 0; y < data.height; ++y) {
			for (int x = 0; x < data.width; ++x) {
				for (int i = 0; i < 10; ++i) {
					//lensSystem.Test(x, y, -1, incomingImageSensor, generatedRay, rand, file);
				}
			}
		}
	}
	return 0;

	double maxSize = 0.5;
	int N = 1;
	for (double y = -maxSize; y <= maxSize; y += maxSize / N) {
		for (double x = -maxSize; x <= maxSize; x += maxSize / N) {
			ray.org = png::vec3(0, 0, -10);
			ray.dir = Normalize(png::vec3(x, y, 1));
			//auto points = Pathtracing(lens, ray);
			//Output(file, points);
		}
	}


	return 0;
}
*/

int main(int argc, char* argv[]) {
	return png::main2D();
	std::string jsonFileTest = "settingData.json";
	auto data = png::LoadData(jsonFileTest);
	data.data.width = 5;
	data.data.height = 5;
	png::LensSystem lensSystem("camera/convex.json", data.data);

	png::Ray ray;
	png::Ray incomingImageSensor;
	std::optional<png::Ray> generatedRay;
	png::Random rand;
	if (0 == 0) {
		std::ofstream file36("testpath.txt");
		for (int i = 0; i < 100; ++i) {
			while (true) {
				lensSystem.Test(0, 2, 0, -1, true, incomingImageSensor, generatedRay, rand, file36);
				if (generatedRay.has_value()) {
					break;
				}
			}
		}

		for (double z = 0; z >= -50; z -= 5.0) {
			std::ofstream file("path" + std::to_string(int(z)) + ".txt");
			for (int i = 0; i < 100; ++i) {
				while (true) {
					lensSystem.Test(2, 2, z, -1, false, incomingImageSensor, generatedRay, rand, file);
					if (generatedRay.has_value()) {
						break;
					}
				}
			}
		}
	}
	else {
		// ray.dir = (0, 0, 1)
		std::optional<png::Ray> generatedRay;
		std::ofstream file("path.txt");
		for (double x = -1.0; x <= 1.0; x += 0.25) {
			lensSystem.RayDir(x, 0.0, -50.0, -1, generatedRay, rand, file);
		}
	}

	double maxSize = 0.5;
	int N = 1;
	for (double y = -maxSize; y <= maxSize; y += maxSize / N) {
		for (double x = -maxSize; x <= maxSize; x += maxSize / N) {
			ray.org = png::vec3(0, 0, -10);
			ray.dir = Normalize(png::vec3(x, y, 1));
			//auto points = Pathtracing(lens, ray);
			//Output(file, points);
		}
	}
	//return 0;

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
	renderer.Render("result");

	return 0;
}
