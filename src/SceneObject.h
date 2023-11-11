#pragma once

#include "Ray.h"
#include "Random.h"
#include "HitRecord.h"
#include "Texture.h"
#include "Transparent.h"

namespace png {
    struct HitRecord;
	struct SettingData;
	double kd(const vec3& color);
	void RayCasting(const Ray& ray, const SettingData& data, int& objectIndex, double& distance);

	enum class PDF_TYPE {
		UniformRandom,
		CosImportSample,
		LightSample,
		Mixture,
	};

	class Material {
	public:
		Material(Texture* color, Texture* emission);
		virtual Ray ScatteredRay(const Ray refRay, HitRecord& rec, const SettingData& data, double& out_pdfWeight, const double spectrum, Random& rand) const = 0;
		vec3 color(double x, double y) const;
		vec3 emission(double x, double y) const;
	public:
		Texture* m_color;
		Texture* m_emission;
	};
	class RefractionMaterial : public Material {
	public:
		RefractionMaterial(const TransparentMaterialType transparentMaterialType, Texture* color, Texture* emission);
		Ray ScatteredRay(const Ray refRay, HitRecord& rec, const SettingData& data, double& out_pdfWeight, const double spectrum, Random& rand) const;
	private:
		TransparentMaterialType m_transparentMaterialType;
	};
	class DiffuseMaterial : public Material {
	public:
		DiffuseMaterial(Texture* color, Texture* emission, PDF_TYPE pdf_type);
		Ray ScatteredRay(const Ray refRay, HitRecord& rec, const SettingData& data, double& out_pdfWeight, const double spectrum, Random& rand) const;
	private:
		PDF_TYPE m_PDF_TYPE;
	};

	class SceneObject {
	public:
		SceneObject(Material* mat);
		bool Hitable(const Ray& ray) const;
		virtual double HitDistance(const Ray& ray) const = 0;
		virtual Ray ScatteredRay(const Ray& refRay, HitRecord& hitrecord, const SettingData& data, double& out_pdfWeight, const double spectrum, Random& rand) const = 0;
		virtual vec3 RandomSurfacePoint(Random& rand) const = 0;
		virtual vec3 Normal(const vec3& point) const = 0;

		//color
		virtual vec3 color(const vec3& point) const = 0;
		virtual vec3 emission(const vec3& point) const = 0;
	public :
		Material* m_material;
	};
	class SphereObject : public SceneObject {
	public:
		SphereObject(vec3 posi, float size, Material* mat);
		double HitDistance(const Ray& ray) const;
		Ray ScatteredRay(const Ray& refRay, HitRecord& hitrecord, const SettingData& data, double& out_pdfWeight, const double spectrum, Random& rand) const;
		vec3 RandomSurfacePoint(Random& rand) const;
		vec3 Normal(const vec3& point) const;

		//color
		vec3 color(const vec3& point) const;
		vec3 emission(const vec3& point) const;

		void SetPosition(vec3 posi);
		void SetSize(float size);

		vec3 position() const;
		float size() const;
		
	private:
		vec3 m_position;
		float m_size;
	};

	class MeshObject : public SceneObject {
	public:
		MeshObject(std::string fileName, Material* mat);
		MeshObject(std::vector<vec3> mesh, Material* mat);
		double HitDistance(const Ray& ray) const;
		Ray ScatteredRay(const Ray& refRay, HitRecord& hitrecord, const SettingData& data, double& out_pdfWeight, const double spectrum, Random& rand) const;
		vec3 RandomSurfacePoint(Random& rand) const;
		vec3 Normal(const vec3& point) const;

		//color
		vec3 color(const vec3& point) const;
		vec3 emission(const vec3& point) const;

	protected:
		std::vector<vec3> m_mesh;
		double m_SumArea;
		std::vector<double> m_area;
	};

	class BoxObject : public MeshObject {
	public :
		BoxObject(const vec3& offset, const vec3& size, Material* mat);
		vec3 RandomSurfacePoint(Random& rand) const;
		vec3 Normal(const vec3& point) const;

		//color
		vec3 color(const vec3& point) const;
		vec3 emission(const vec3& point) const;
	private:
	};

	class PlaneObject : public MeshObject {
	public:
		PlaneObject(const vec3& offset, const vec3& up, const vec3& right, Material* mat);
		vec3 RandomSurfacePoint(Random& rand) const;
		vec3 Normal(const vec3& point) const;

		//color
		vec3 color(const vec3& point) const;
		vec3 emission(const vec3& point) const;
	private:
		vec3 m_up, m_right;
		vec3 m_offset;
	};
}