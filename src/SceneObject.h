#pragma once

#include "Ray.h"
#include "Random.h"
#include "HitRecord.h"
#include "Texture.h"
#include "Transparent.h"

namespace png {
    struct HitRecord;
	double kd(const vec3& color);

	class Material {
	public:
		Material(Texture* color, Texture* emission);
		virtual Ray ScatteredRay(const Ray refRay, HitRecord& rec, const double spectrum, Random& rand) const = 0;
		vec3 color(double x, double y) const;
		vec3 emission(double x, double y) const;
	private:
		Texture* m_color;
		Texture* m_emission;
	};
	class RefractionMaterial : public Material {
	public:
		RefractionMaterial(const TransparentMaterialType transparentMaterialType, Texture* color, Texture* emission);
		Ray ScatteredRay(const Ray refRay, HitRecord& rec, const double spectrum, Random& rand) const;
	private:
		TransparentMaterialType m_transparentMaterialType;
	};
	class DiffuseMaterial : public Material {
	public:
		DiffuseMaterial(Texture* color, Texture* emission);
		Ray ScatteredRay(const Ray refRay, HitRecord& rec, const double spectrum, Random& rand) const;
	};

	class SceneObject {
	public:
		SceneObject(Material* mat);
		bool Hitable(const Ray& ray) const;
		virtual double HitDistance(const Ray& ray) const = 0;
		virtual Ray ScatteredRay(const Ray& refRay, HitRecord& hitrecord, const double spectrum, Random& rand) const = 0;

		//color
		virtual vec3 color(const vec3& point) const = 0;
		virtual vec3 emission(const vec3& point) const = 0;
	protected :
		Material* material;
	};
	class SphereObject : public SceneObject {
	public:
		SphereObject(vec3 posi, float size, Material* mat);
		double HitDistance(const Ray& ray) const;
		Ray ScatteredRay(const Ray& refRay, HitRecord& hitrecord, const double spectrum, Random& rand) const;

		//color
		vec3 color(const vec3& point) const;
		vec3 emission(const vec3& point) const;

		void SetPosition(vec3 posi);
		void SetSize(float size);

		vec3 position() const;
		float size() const;
		
	protected:
		vec3 m_position;
		float m_size;
	};

	class SpotlightObject : public SphereObject {
	public:
		SpotlightObject(vec3 posi, float size, vec3 front, double dotVal, Material* mat);

		double HitDistance(const Ray& ray) const override;

	private:
		vec3 m_front;
		double m_dotVal;
	};

	class MeshObject : public SceneObject {
	public:
		MeshObject(std::string fileName, Material* mat);
		MeshObject(std::vector<vec3> mesh, Material* mat);
		double HitDistance(const Ray& ray) const;
		Ray ScatteredRay(const Ray& refRay, HitRecord& hitrecord, const double spectrum, Random& rand) const;

		//color
		vec3 color(const vec3& point) const;
		vec3 emission(const vec3& point) const;

	protected:
		std::vector<vec3> mesh;
	};

	class BoxObject : public MeshObject {
	public :
		BoxObject(const vec3& offset, const vec3 size, Material* mat);

		//color
		vec3 color(const vec3& point) const;
		vec3 emission(const vec3& point) const;
	private:
	};

	class PlaneObject : public MeshObject {
	public:
		PlaneObject(const vec3& offset, const vec3& up, const vec3& right, Material* mat);

		//color
		vec3 color(const vec3& point) const;
		vec3 emission(const vec3& point) const;
	private:
		vec3 m_up, m_right;
		vec3 m_offset;
	};
}