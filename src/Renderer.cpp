#include "Renderer.h"

#include <numbers>
#include <random>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb_image_write.h>
#include <iostream>

namespace png {
	namespace {
		constexpr double PI = 3.14159265358979323846;
		template <typename T>
		struct nullable {
			T value;
			bool has_value = false;
		};
	}
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

	nullable<double> intersect(const Ray& ray, const Object& obj) {
		const vec3 p_o = obj.position - ray.org;
		const double b = Dot(p_o, ray.dir);
		const double D4 = b * b - Dot(p_o, p_o) + obj.size * obj.size;

		if (D4 < 0.0)
			return nullable<double>{0, false};

		const double sqrt_D4 = sqrt(D4);
		const double t1 = b - sqrt_D4, t2 = b + sqrt_D4;

		const float minValue = 1e-5;
		if (t1 < minValue && t2 < minValue)
			return nullable<double>{0, false};

		if (t1 > 0.001) {
			return nullable<double>{t1, true};
		}
		else {
			return nullable<double>{t2, true};
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

	vec3 PathTracing(const Ray ray, SettingData& data, std::random_device& rand) {
		int hitObject = -1;
		double dis = std::numeric_limits<double>::max();
		{
			for (int i = 0; i < data.object.size(); ++i) {
				auto& obj = data.object[i];
				auto tmp_dis = intersect(ray, obj);
				//float tmp_dis = hit_sphere(obj.position, obj.size, ray);
				if (tmp_dis.has_value && tmp_dis.value < dis && tmp_dis.value > 0) {
					dis = tmp_dis.value;
					hitObject = i;
				}
			}
		}

		if (hitObject != -1) {
			auto& obj = data.object[hitObject];
			if (dis > 0) {
				if (random(rand) <= obj.material->kd()) {
					const auto hitPoint = ray.dir * dis + ray.org;
					const auto normal_hitedPoint = Normalize(hitPoint - obj.position);
					auto nextRay = obj.material->ScatteredRay();
					/*
					const auto orienting_normal = Dot(normal_hitedPoint, ray.dir) < 0.0
						? normal_hitedPoint : (normal_hitedPoint * -1.0);
					nextRay.org = hitPoint;
					vec3 u, v, w;
					w = orienting_normal;
					const auto r1 = 2 * PI * random(rand);
					const auto r2 = random(rand);
					const auto r2s = sqrt(r2);
					if (fabs(w.x) > std::numeric_limits<float>::min()) {
						u = Normalize(Cross(vec3(0, 1, 0), w));
					}
					else {
						u = Normalize(Cross(vec3(1, 0, 0), w));
					}
					v = Cross(w, u);
					nextRay.dir = Normalize((
						u * cos(r1) * r2s +
						v * sin(r1) * r2s +
						w * sqrt(1.0 - r2))
					);
					*/
					auto nextPathTracing = PathTracing(nextRay, data, rand);
					return obj.material->colorKD() * nextPathTracing + obj.material->emission;
				}
				else {
					return obj.material->emission;
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
				vec3 accumulatedColor = vec3(0,0,0);
				for (int sx = 1; sx <= data.superSamples; ++sx) {
					for (int sy = 1; sy <= data.superSamples; ++sy) {
						for (int s = 0; s < data.samples; ++s) {
							const float rate = 1.0 / (1 + data.superSamples);
							vec3 dir = Normalize(
								l_camX * fovx * (2.0f * ((double)x + rate * sx) / data.width - 1.0f) +
								l_camY * fovy * (2.0f * ((double)y + rate * sy) / data.height - 1.0f) +
								l_camZ
							);
							auto cal = PathTracing(Ray(data.camera.origin, dir), data, rnd) / data.superSamples / data.superSamples / data.samples;
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