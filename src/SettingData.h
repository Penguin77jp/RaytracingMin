#pragma once
#include <string>
#include <cmath>
#include <nlohmann/json.hpp>
#include "Ray.h"

namespace png {
	struct Material {
	public:
		vec3 color;
		vec3 emission;
		vec3 colorKD() const {
			return color / kd();
		}
		float kd() const {
			return std::max(std::max(color.x, color.y), color.z);
		}
	};
	struct RefractionMaterial : public Material {
	};
	struct DiffuseMaterial : public Material {
	};

	class Object {
	public:
		Material* material;
		bool Hitable(const Ray& ray) const { return HitDistance(ray) > 0; }
		virtual double HitDistance(const Ray& ray) const = 0;
		virtual Ray ScatteredRay() const = 0;
	};
	class SphereObject : public Object {
	public:
		double HitDistance(const Ray& ray) const {
			const vec3 p_o = m_position - ray.org;
			const double b = Dot(p_o, ray.dir);
			const double D4 = b * b - Dot(p_o, p_o) + m_size * m_size;

			if (D4 < 0.0)
				return 0;

			const double sqrt_D4 = sqrt(D4);
			const double t1 = b - sqrt_D4, t2 = b + sqrt_D4;

			const float minValue = 1e-5;
			if (t1 < minValue && t2 < minValue)
				return 0;

			if (t1 > 0.001) {
				return t1;
			}
			else {
				return t2;
			}
		}
	private:
		vec3 m_position;
		float m_size;
	};

	struct Camera {
		vec3 origin, target, upVec;
		float fov;
	};
	struct SettingData {
		int width, height, samples, superSamples;
		Camera camera;
		std::vector<Object*> object;
	};

	struct LoadData {
		LoadData(std::string);
		static void SaveSampleJson(std::string);

		static void to_json(nlohmann::json& json, const SettingData& data);
		static void from_json(const nlohmann::json& json, SettingData& data);

		SettingData data;
	};
}