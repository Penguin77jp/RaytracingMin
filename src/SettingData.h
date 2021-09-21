#pragma once
#include <string>
#include <cmath>
#include <nlohmann/json.hpp>
#include "Ray.h"

namespace png {
	class Material {
	public:
		Material() = default;
		vec3 color;
		vec3 emission;
		vec3 colorKD() const {
			return color / kd();
		}
		float kd() const {
			return std::max(std::max(color.x, color.y), color.z);
		}
		virtual Ray ScatteredRay
		(
			const Ray& inScatteredRay,
			std::random_device& randomDevice
		)const = 0;
	};
	class RefractionMaterial : public Material {
	public:
		RefractionMaterial() = default;
		Ray ScatteredRay
		(
			const Ray& inScatteredRay,
			std::random_device& randomDevice
		)const {
			return Ray();
		}
	};
	class DiffuseMaterial : public Material {
	public:
		DiffuseMaterial()
			:Material()
		{}
		Ray ScatteredRay
		(
			const Ray& inScatteredRay,
			std::random_device& randomDevice
		)const {
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
		}
	};

	struct Object {
		vec3 position;
		float size;
		Material* material;
	};
	struct Camera {
		vec3 origin, target, upVec;
		float fov;
	};
	struct SettingData {
		int width, height, samples, superSamples;
		Camera camera;
		std::vector<Object> object;
	};

	struct LoadData {
		LoadData(std::string);
		static void SaveSampleJson(std::string);

		static void to_json(nlohmann::json& json, const SettingData& data);
		static void from_json(const nlohmann::json& json, SettingData& data);

		SettingData data;
	};
}