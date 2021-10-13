#include "SceneObject.h"
#include "Constant.h"


namespace png {
	Material::Material(vec3 color, vec3 emission)
		: m_color(color)
		, m_emission(emission)
	{}
	vec3 Material::colorKD() const {
		return color() / kd();
	}
	float Material::kd() const {
		return std::max(std::max(color().x, color().y), color().z);
	}
	vec3 Material::color() const { return m_color; }
	vec3 Material::emission() const { return m_emission; }

	RefractionMaterial::RefractionMaterial(vec3 color, vec3 emission)
		: Material(color, emission)
	{}
	Ray RefractionMaterial::ScatteredRay(const Ray refRay, const vec3 hitPoint, const vec3 hitedNormal, const double spectrum, std::random_device& rand) const {
		return Ray();
	}

	DiffuseMaterial::DiffuseMaterial(vec3 color, vec3 emission)
		: Material(color, emission)
	{}
	Ray DiffuseMaterial::ScatteredRay(const Ray refRay, const vec3 hitPoint, const vec3 hitedNormal, const double spectrum, std::random_device& rand) const {
		Ray nextRay;
		return nextRay;
		/*
		const auto normalHitedPoint = Normalize(hitPoint - this->position());
		const auto orienting_normal = Dot(normalHitedPoint, refRay.dir) < 0.0
			? normalHitedPoint : (normalHitedPoint * -1.0);
		nextRay.org = hitPoint;
		vec3 u, v, w;
		w = orienting_normal;
		const auto r1 = 2 * PI * RandomGenerate(rand);
		const auto r2 = RandomGenerate(rand);
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
		return nextRay;
		*/
	}

	SceneObject::SceneObject(Material* mat)
	: material(mat)
	{}


	SphereObject::SphereObject(vec3 posi, float size, Material* mat)
	: m_position(posi),
		m_size(size),
		SceneObject(mat)
	{}

	bool SceneObject::Hitable(const Ray& ray) const { return HitDistance(ray) > 0; }
	double SphereObject::HitDistance(const Ray& ray) const {
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
	Ray SphereObject::ScatteredRay(const Ray refRay, const double spectrum, std::random_device& rand) const {
		const auto distance = HitDistance(refRay);
		const auto hitPoint = refRay.dir * distance + refRay.org;
		return material->ScatteredRay(refRay, hitPoint, spectrum, rand);
		return Ray();
	}
	//vec3 SphereObject::Normal(vec3) const { return vec3(); }

	void SphereObject::SetPosition(vec3 posi) { m_position = posi; }
	void SphereObject::SetSize(float size) { m_size = size; }

	vec3 SphereObject::position() const { return m_position; }
	float SphereObject::size() const { return m_size; }

}