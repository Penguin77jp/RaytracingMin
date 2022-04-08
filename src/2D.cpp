#pragma once

#include <optional>
#include <fstream>

#include "2D.h"
#include "Camera.h"
#include "Ray.h"
#include "Random.h"
#include "stb_image_write.h"

namespace png {
	std::optional<png::vec3> Refract2D(png::vec3 incomingDir, png::vec3 lensNormal, double n1, double n2) {
		auto cos2 = 1.0 - pow2(n1 / n2) * (1.0 - pow2(Dot(-incomingDir, lensNormal)));
		if (cos2 <= 0.0) {
			// full refract
			return std::nullopt;
		}
		return n1 / n2 * incomingDir - (n1 / n2 * Dot(incomingDir, lensNormal) + std::sqrt(cos2)) * lensNormal;
	}

	void Intersect2D(png::vec3 sphereCenter, LensSurface& lensData, png::Ray ray, std::optional<png::vec3>& hitpoint, png::vec3& normal) {
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
			if (pow2(hitpoint.value().x) + pow2(hitpoint.value().y) > pow2(lensData.h)) {
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
			hitpoint = ray.org + t * ray.dir;
			if (pow2(hitpoint.value().x) + pow2(hitpoint.value().y) > pow2(lensData.h)) {
				hitpoint = std::nullopt;
				return;
			}

			normal = png::vec3(0, 0, -1);
			if (Dot(ray.dir, normal) > 0) {
				normal = -normal;
			}
		}
	}

	void GenerateRay(int x, int y, int superSampleX, int superSampleY, double spectrum, Ray& rayIncomingSensor, std::optional<Ray>& generatedRay, Random& rand, const LensSystem& lens, std::vector<vec3>& points, const SettingData& data) {
		Ray lensRay;
		double z = 0.0;
		{
			// assumption thin lens
			auto yover = 1.0 / lens.m_lensSurface[0].f - 1.0 / lens.m_targetDistance;
			z = 1.0 / yover;

			/*
			* assumption thick lens
			auto Si = -std::min(m_lensSurface[1].pp, 0.0) + m_targetDistance;
			auto So = m_lensSurface[0].f * Si / (Si - m_lensSurface[0].f);
			z = So - std::max(m_lensSurface[0].pp, 0.0);
			*/
		}
		{
			const auto rate = 1.0 / (1 + lens.superSamples());
			const auto maxImageSensorResolution = std::max(lens.m_data.width, lens.m_data.height);
			const auto indexX = (double)x;
			const auto indexY = (double)y;
			lensRay.org = vec3(0.0
				, (2.0 * indexY / lens.m_data.height - 1.0) * lens.m_data.height
				, lens.m_lensSurface[0].z - z);
			std::optional<Ray> outRay;
			//double xi_theta = std::acos(std::sqrt(rand.RandomGenerate()));
			//double xi_phi = 2.0 * png::PI * rand.RandomGenerate();
			//lensRay.dir = Normalize(vec3(std::sin(xi_theta) * std::cos(xi_phi), std::sin(xi_theta) * std::sin(xi_phi), std::cos(xi_theta)));
			lensRay.dir = Normalize(vec3(0, 2.0*rand.RandomGenerate()-1.0, 1.0));

			rayIncomingSensor = lensRay;
		}

		points.push_back(lensRay.org + data.cameraOrigin);
		{
			double n1 = 1.0;
			for (int i = 0; i < lens.m_lensSurface.size(); ++i) {
				auto targetLens = lens.m_lensSurface[i];

				std::optional<png::vec3> hitpoint;
				png::vec3 lensNormal;
				Intersect2D(png::vec3(0, 0, targetLens.z + targetLens.r), targetLens, lensRay, hitpoint, lensNormal);
				if (!hitpoint.has_value()) {
					// no ray
					generatedRay = std::nullopt;
					return;
				}

				auto refractedDir = Refract2D(lensRay.dir, lensNormal, n1, targetLens.ior);

				if (!refractedDir.has_value()) {
					// no ray
					generatedRay = std::nullopt;
					return;
				}
				lensRay.dir = refractedDir.value();
				lensRay.org = hitpoint.value();
				points.push_back(lensRay.org + data.cameraOrigin);
				n1 = targetLens.ior;
			}
		}

		auto m_camX = vec3(1, 0, 0);
		auto m_camY = vec3(0, 1, 0);
		auto m_camZ = vec3(0, 0, 1);
		auto raySize = 1.0;
		Ray generatedRayWorldSpace;
		generatedRayWorldSpace.org = m_camX * raySize * lensRay.org.x + m_camY * raySize * lensRay.org.y + data.cameraOrigin;
		generatedRayWorldSpace.dir = Normalize(m_camX * raySize * lensRay.dir.x + m_camY * raySize * lensRay.dir.y + m_camZ * lensRay.dir.z);
		generatedRay = generatedRayWorldSpace;
	}

	vec3 PathTrace(const int y, Random& rand, LoadData& loadedData, LensSystem& lens, std::vector<vec3>& points) {
		std::optional<Ray> ray;

		png::Ray incomingImageSensor;
		std::vector<vec3> tmp_points;
		while (true) {
			GenerateRay(0, y, 0, 0, -1, incomingImageSensor, ray, rand, lens, tmp_points, loadedData.data);
			if (ray.has_value()) {
				break;
			}
			tmp_points.clear();
		}
		for (int i = 0; i < tmp_points.size(); ++i) {
			points.push_back(tmp_points[i]);
		}

		if (ray.has_value()) {
			auto focalZ = loadedData.data.cameraTarget.z;
			auto t = (focalZ - ray.value().org.z) / ray.value().dir.z;
			auto point = ray.value().org + t * ray.value().dir;
			points.push_back(point);

			auto objectSize = 0.1;
			if (point.x > objectSize) {
				return vec3(1, 0, 0);
			}
			else if (point.x > -objectSize) {
				return vec3(0, 1, 0);
			}
			else {
				return vec3(0, 0, 1);
			}
		}

		return vec3();
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

	int main2D() {

		LoadData loadedData("settingData.json");
		loadedData.data.height = 3;
		loadedData.data.samples = 100;
		const auto height = loadedData.data.height;
		const auto sample = loadedData.data.samples;
		LensSystem lens("camera/convex.json", loadedData.data);

		Random rand;
		std::vector<unsigned char> image(height * 3);
		std::vector<vec3> points;
		for (int i = 0; i < 2; ++i) {
			PathTrace(0, rand, loadedData, lens, points);
		}
		std::ofstream output("2d.txt");
		Output(output, points);
		return 0;
		for (int y = 0; y < height; ++y) {
			auto color = vec3();
			for (int s = 0; s < sample; ++s) {
				color += PathTrace(y, rand, loadedData, lens, points) / sample;
			}
			image[y * 3 + 0] = 255 * color.x;
			image[y * 3 + 1] = 255 * color.y;
			image[y * 3 + 2] = 255 * color.z;
		}

		// output
		const auto size = 50;
		std::vector<unsigned char> sizedImage(size * height * size * 3);
		for (int y = 0; y < size * height; ++y) {
			for (int x = 0; x < size; ++x) {
				auto index = (y / size) * 3 + (x / size) * 3;
				sizedImage[y * size * 3 + x * 3 + 0] = image[index + 0];
				sizedImage[y * size * 3 + x * 3 + 1] = image[index + 1];
				sizedImage[y * size * 3 + x * 3 + 2] = image[index + 2];
			}
		}
		stbi_write_bmp("2d.bmp", size, height * size, 3, sizedImage.data());

		return 0;
	}
}