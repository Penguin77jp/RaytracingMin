#pragma once
#include <string>
#include <cmath>
#include <nlohmann/json.hpp>
#include "Ray.h"

namespace png {
	struct Material {
		vec3 color;
		vec3 emission;
		vec3 colorKD() const {
			return color / kd();
		}
		float kd() const {
			return std::max(std::max(color.x, color.y), color.z);
		}
	};
	struct Object {
		Object(const Material mat);
		virtual vec3 ComputeSurfacePoint(const std::function<double()> randGene) const = 0;
		virtual bool Intersect(const Ray& ray, double& out_dis, vec3 out_normal) const = 0;

		Material m_material;
	};
	struct SphereObject : public Object {
		SphereObject(const vec3 position, const float size, const Material mat);
		vec3 ComputeSurfacePoint(const std::function<double()> randGene) const;
		bool Intersect(const Ray& ray, double& out_dis, vec3 out_normal) const;

		vec3 m_position;
		float m_size;
	};
	struct PlaneObject : public Object{
		PlaneObject(const vec3 position, const vec3 up, const vec3 target, const double width, const Material mat);
		vec3 ComputeSurfacePoint(const std::function<double()> randGene) const;
		bool Intersect(const Ray& ray, double& out_dis, vec3 out_normal) const;

		vec3 m_position;
		vec3 m_up;
		vec3 m_right;
		vec3 m_normal;
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