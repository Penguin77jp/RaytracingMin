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

	double random(std::mt19937& gene) {
		return (double)gene() / std::numeric_limits<unsigned int>::max();
	}

	vec3 PathTracing(const Ray ray, SettingData& data, RandType& rand) {
		for (int i = 0; i < data.object.size(); ++i) {
			auto& obj = data.object[i];
			float dis = hit_sphere(obj.position, obj.size, ray);
			if (dis > 0) {
				if (random(rand) <= obj.material.kd) {
					auto tmpRay = ray;
					tmpRay.org = ray.dir * dis + ray.org;
					auto tmpDir = Normalize(tmpRay.org - obj.position);
					while (true) {
						auto randVec3 = Normalize(vec3(
							2.0 * random(rand) - 1.0,
							2.0 * random(rand) - 1.0,
							2.0 * random(rand) - 1.0
						));
						if (Dot(tmpDir, randVec3) >= 0) {
							tmpDir = randVec3;
							break;
						}
					}
					tmpRay.dir = tmpDir;

					auto nextPathTracing = PathTracing(tmpRay, data, rand);
					auto dot = -Dot(ray.dir, Normalize(tmpRay.org - obj.position));
					return obj.material.color * nextPathTracing * dot / obj.material.kd + obj.material.emission;
				}
				else {
					return obj.material.emission;
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
		double fovx = std::cos((90.0f - data.camera.fov / 2.0f) / 180 * std::numbers::pi);
		double fovy = fovx * data.height / data.width;
		//random
		std::random_device rnd;
		RandType mt(rnd());

		for (int y = 0; y < data.height; ++y) {
			std::cout << y << " / " << data.height << std::endl;
#ifdef _DEBUG
#else
#pragma omp parallel for
#endif

			for (int x = 0; x < data.width; ++x) {
				for (int s = 0; s < data.sample; ++s) {
					vec3 dir = Normalize(
						l_camX * fovx * (2.0f * ((float)x + random(mt)) / data.width - 1.0f) +
						l_camY * fovy * (2.0f * ((float)y + random(mt)) / data.height - 1.0f) +
						l_camZ
					);
					auto cal = PathTracing(Ray(data.camera.origin, dir), data, mt) / data.sample;
					image[x * 3 + y * data.width * 3] += cal.x;
					image[x * 3 + y * data.width * 3 + 1] += cal.y;
					image[x * 3 + y * data.width * 3 + 2] += cal.z;
				}
			}
		}
		auto resultImage = std::vector<unsigned char>(data.width * data.height * 3);
		for (int i = 0; i < resultImage.size(); ++i) {
			resultImage[i] = (unsigned char)255 * std::min(image[i],1.0);
		}

		stbi_write_bmp((fileName + ".bmp").c_str(), data.width, data.height, 3, resultImage.data());
	}
}