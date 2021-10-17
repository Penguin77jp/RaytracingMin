#include "Renderer.h"

#include <numbers>
#include <random>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb_image_write.h>
#include <iostream>

namespace png {
	using RandType = std::mt19937;
	Renderer::Renderer(SettingData& data)
		:image(std::vector<double>(data.width* data.height * 3))
		, data(data)
	{}

	float hit_sphere(const vec3 center, double radius, const Ray r) {
		vec3 oc = r.org - center;
		auto a = Dot(r.dir, r.dir);
		auto b = 2.0 * Dot(oc, r.dir);
		auto c = Dot(oc, oc) - radius * radius;
		auto discriminant = b * b - 4 * a * c;
		if (discriminant < 0) {
			return -1.0;
		}
		else {
			return (-b - sqrt(discriminant)) / (2.0 * a);
		}
	}

	double random(std::random_device& gene) {
		return (double)gene() / std::numeric_limits<unsigned int>::max();
	}

	vec3 clampColor(const vec3& ref, double min, float max) {
		return vec3(
			std::min<double>(std::max<double>(ref.x, min), max),
			std::min<double>(std::max<double>(ref.y, min), max),
			std::min<double>(std::max<double>(ref.z, min), max)
		);
	}

	void PrintfDebug(std::string log) {
#ifdef _DEBUG
		std::cout << log;
#endif
	}

	vec3 PathTracing(const Ray ray, const double spectrum, SettingData& data, std::random_device& rand) {
		int hitObject = -1;
		double dis = std::numeric_limits<double>::max();
		{
			for (int i = 0; i < data.object.size(); ++i) {
				auto& obj = data.object[i];
				auto tmp_dis = obj->HitDistance(ray);
				//float tmp_dis = hit_sphere(obj.position, obj.size, ray);
				if (tmp_dis < dis && tmp_dis > 0) {
					dis = tmp_dis;
					hitObject = i;
				}
			}
		}

		if (hitObject != -1) {
			auto& obj = data.object[hitObject];
			if (dis > 0) {
				if (random(rand) <= obj->material->kd()) {
					auto nextRay = obj->ScatteredRay(ray, 0, rand);
					auto nextPathTracing = PathTracing(nextRay,spectrum, data, rand);
					return obj->material->color() / obj->material->kd() * nextPathTracing + obj->material->emission();
				}
				else {
					return obj->material->emission();
				}
			}
		}
		return vec3();
	}

	void Renderer::Render(std::string fileName) {
		//dir
		auto direction = Normalize(data.camera.target - data.camera.origin);
		auto l_camX = -Normalize(Cross(direction, data.camera.upVec));
		auto l_camY = Cross(l_camX, direction);
		auto l_camZ = direction;
		//fov
		double fovx = data.camera.fov;
		double fovy = fovx * data.height / data.width;
		//random
		std::random_device rnd;

		for (int y = 0; y < data.height; ++y) {
			std::cout << y << " / " << data.height << std::endl;
#ifdef _DEBUG
#else
#pragma omp parallel for
#endif

			for (int x = 0; x < data.width; ++x) {
				vec3 accumulatedColor = vec3(0, 0, 0);
				for (int sx = 1; sx <= data.superSamples; ++sx) {
					for (int sy = 1; sy <= data.superSamples; ++sy) {
						for (int s = 0; s < data.samples; ++s) {
							const float rate = 1.0 / (1 + data.superSamples);
							vec3 dir = Normalize(
								l_camX * fovx * (2.0f * ((double)x + rate * sx) / data.width - 1.0f) +
								l_camY * fovy * (2.0f * ((double)y + rate * sy) / data.height - 1.0f) +
								l_camZ
							);
							const double spectrum = random(rnd)*440.0 + 390;
							auto cal = PathTracing(Ray(data.camera.origin, dir), spectrum, data, rnd) / data.superSamples / data.superSamples / data.samples;
							cal = clampColor(cal, 0, 1);
							accumulatedColor += cal;
						}
					}
				}
				image[x * 3 + y * data.width * 3] += accumulatedColor.x;
				image[x * 3 + y * data.width * 3 + 1] += accumulatedColor.y;
				image[x * 3 + y * data.width * 3 + 2] += accumulatedColor.z;
			}
		}
		auto resultImage = std::vector<unsigned char>(data.width * data.height * 3);
		for (int i = 0; i < resultImage.size(); ++i) {
			resultImage[i] = (unsigned char)255 * std::min(image[i], 1.0);
		}

		stbi_write_jpg((fileName + ".jpg").c_str(), data.width, data.height, 3, resultImage.data(), 60);
		stbi_write_bmp((fileName + ".bmp").c_str(), data.width, data.height, 3, resultImage.data());
	}
}