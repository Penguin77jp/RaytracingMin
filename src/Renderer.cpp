#include "Renderer.h"

#include <numbers>
#include <random>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb_image_write.h>
#include <iostream>
#include <numbers>

namespace png {
	namespace {
		constexpr double PI = std::numbers::pi;
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

	vec3 PathTracing(const Ray ray, SettingData& data, const std::function<double()>& randomGene) {
		int hitObject = -1;
		vec3 normal;
		double dis = std::numeric_limits<double>::max();
		{
			for (int i = 0; i < data.object.size(); ++i) {
				auto obj = data.object[i];
				double tmp_dis;
				vec3 tmp_normal;
				auto tmp_intersect = obj->Intersect(ray, tmp_dis, tmp_normal);
				if (tmp_intersect && tmp_dis < dis && tmp_dis > 0) {
					dis = tmp_dis;
					normal = tmp_normal;
					hitObject = i;
				}
			}
		}

		if (hitObject != -1) {
			auto* obj = data.object[hitObject];
			if (dis > 0) {
				if (randomGene() <= obj->m_material.kd()) {
					const auto hitPoint = ray.dir * dis + ray.org;
					const auto normal_hitedPoint = normal;
					//const auto normal_hitedPoint = Normalize(hitPoint - obj.position);
					const auto orienting_normal = Dot(normal_hitedPoint, ray.dir) < 0.0
						? normal_hitedPoint : (normal_hitedPoint * -1.0);
					auto nextRay = ray;
					nextRay.org = hitPoint;
					vec3 u, v, w;
					w = orienting_normal;
					const auto r1 = 2 * PI * randomGene();
					const auto r2 = randomGene();
					const auto r2s = sqrt(r2);
					/*
					error
					const auto randPhi = 1.0 / 2 / std::numbers::pi * random(rand);
					const auto randTheta = std::asin(std::sqrt(random(rand)));
					*/
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
					auto nextPathTracing = PathTracing(nextRay, data, rand);
					return obj->m_material.colorKD() * nextPathTracing + obj->m_material.emission;
				}
				else {
					return obj->m_material.emission;
				}
			}
		}
		return vec3();
	}

	vec3 SurfaeSample(const Ray ray, SettingData& data, const int depth, const std::function<double()>& randomGene) {
		if (depth <= 0) {
			// path check
			vec3 hitpoint;
			vec3 normal;
			int hitObjectIndex = -1;
			double dis = std::numeric_limits<double>::max();
			for (int i = 0; i < data.object.size(); ++i) {
				auto& obj = data.object[i];
				double tmp_dis;
				vec3 tmp_normal;
				auto tmp_intersect = obj->Intersect(ray, tmp_dis, tmp_normal);
				if (tmp_intersect && tmp_dis < dis && tmp_dis > 0) {
					dis = tmp_dis;
					hitpoint = ray.org + ray.dir * dis;
					hitObjectIndex = i;
					normal = tmp_normal;
				}
			}
			if (hitObjectIndex == -1) {
				return vec3();
			}
			if (hitObjectIndex == 1) {
				int hoge = 0;
				hoge = 1;
			}
			const auto& hitObject = data.object[hitObjectIndex];
			auto cal = SurfaeSample(Ray(hitpoint, normal), data, depth + 1, randomGene);
			return cal * hitObject->m_material.color + hitObject->m_material.emission;
		}
		else {
			const int objectIndex = randomGene() * data.object.size();
			const auto object = data.object[objectIndex];
			const auto surfacePoint = object->ComputeSurfacePoint(randomGene);

			if (objectIndex == 0) {
				int hoge = 0;
				hoge = 1;
			}

			// path check
			auto pathCheckRay = Ray(ray.org, Normalize(surfacePoint - ray.org));
			vec3 hitpoint, normal;
			int hitObject = -1;
			{
				double dis = std::numeric_limits<double>::max();
				for (int i = 0; i < data.object.size(); ++i) {
					auto obj = data.object[i];
					double tmp_dis;
					vec3 tmp_normal;
					auto tmp_intersect = obj->Intersect(pathCheckRay, tmp_dis, tmp_normal);
					if (tmp_intersect && tmp_dis < dis && tmp_dis> 0) {
						dis = tmp_dis;
						normal = tmp_normal;
						hitpoint = pathCheckRay.org + pathCheckRay.dir * dis;
						hitObject = i;
					}
				}
				if (hitObject == -1 || hitObject != objectIndex) {
					return vec3();
				}
			}
			const auto dir = Normalize(surfacePoint - ray.org);
			auto dot1 = Dot(ray.dir, dir);
			auto dot2 = Dot(-dir, normal);
			if (dot1 <= 0 || dot2 <= 0) {
				return vec3();
			}
			auto surfacePointDisntace = Magnitude(ray.org - hitpoint);
			auto surfacePointProbabiliy = std::max(surfacePointDisntace, 1.0);
			auto surfacePointDiv = std::min(surfacePointDisntace, 1.0);
			if (dot1 * dot2 * object->m_material.kd() * surfacePointProbabiliy < randomGene()) {
				return object->m_material.emission;
			}
			auto colNext = SurfaeSample(Ray(hitpoint, normal), data, depth + 1, randomGene);
			return colNext * object->m_material.colorKD() * surfacePointDiv + object->m_material.emission;
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
		std::random_device seed;
		std::mt19937 rnd(seed());
		std::uniform_real_distribution<double> dist(0.0, 1.0);
		auto randomGene = [&rnd, &dist]() {return dist(rnd); };

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
							vec3 cal = { 0,0,0 };
							if (randomGene() >= 1.0) {
								cal = SurfaeSample(Ray(data.camera.origin, dir), data, 0, randomGene);
							}
							else {
								cal = PathTracing(Ray(data.camera.origin, dir), data, randomGene);
							}
							cal = cal / data.superSamples / data.superSamples / data.samples;
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