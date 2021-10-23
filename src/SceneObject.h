#pragma once

#include "Ray.h"
#include "Random.h"
#include "HitRecord.h"

namespace png {
	class Material {
	public:
		Material(vec3 color, vec3 emission);
		virtual Ray ScatteredRay(const Ray refRay, HitRecord& rec, const double spectrum, std::random_device& rand) const = 0;
		float kd() const;
		vec3 color() const;
		vec3 emission() const;
	private:
		vec3 m_color;
		vec3 m_emission;
	};
	class RefractionMaterial : public Material {
	public:
		RefractionMaterial(vec3 color, vec3 emission);
		Ray ScatteredRay(const Ray refRay, HitRecord& rec, const double spectrum, std::random_device& rand) const;
	};
	class DiffuseMaterial : public Material {
	public:
		DiffuseMaterial(vec3 color, vec3 emission);
		Ray ScatteredRay(const Ray refRay, HitRecord& rec, const double spectrum, std::random_device& rand) const;
	};

	class SceneObject {
	public:
		SceneObject(Material* mat);
		Material* material;
		bool Hitable(const Ray& ray) const;
		virtual double HitDistance(const Ray& ray) const = 0;
		virtual Ray ScatteredRay(const Ray& refRay, HitRecord& hitrecord, const double spectrum, std::random_device& rand) const = 0;
		//virtual vec3 Normal(vec3 posi) const = 0;
	};
	class SphereObject : public SceneObject {
	public:
		SphereObject(vec3 posi, float size, Material* mat);
		double HitDistance(const Ray& ray) const;
		Ray ScatteredRay(const Ray& refRay, HitRecord& hitrecord, const double spectrum, std::random_device& rand) const;
		//vec3 Normal(vec3) const;

		void SetPosition(vec3 posi);
		void SetSize(float size);

		vec3 position() const;
		float size() const;
		
	private:
		vec3 m_position;
		float m_size;
	};
}