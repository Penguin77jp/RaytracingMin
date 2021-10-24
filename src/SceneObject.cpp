#include "SceneObject.h"
#include "Constant.h"
#include "Color.h"


namespace png {
	Material::Material(vec3 color, vec3 emission)
		: m_color(color)
		, m_emission(emission)
	{}
	float Material::kd() const {
		return std::max(std::max(color().x, color().y), color().z);
	}
	vec3 Material::color() const { return m_color; }
	vec3 Material::emission() const { return m_emission; }

	RefractionMaterial::RefractionMaterial(vec3 color, vec3 emission)
		: Material(color, emission)
	{}

	static double reflectance(double cosine, double ref_idx) {
		// Use Schlick's approximation for reflectance.
		auto r0 = (1 - ref_idx) / (1 + ref_idx);
		r0 = r0 * r0;
		return r0 + (1 - r0) * pow((1 - cosine), 5);
	}

	inline vec3 reflect(const vec3& v, const vec3& n) {
		return v - 2 * Dot(v, n) * n;
	}

	double length_squared(const vec3& e) {
		return e.x * e.x + e.y * e.y + e.z * e.z;
	}

	inline vec3 refract(const vec3& uv, const vec3& n, double etai_over_etat) {
		auto cos_theta = fmin(Dot(-uv, n), 1.0);
		vec3 r_out_perp = etai_over_etat * (uv + cos_theta * n);
		vec3 r_out_parallel = -sqrt(fabs(1.0 - length_squared(r_out_perp))) * n;
		return r_out_perp + r_out_parallel;
	}

	Ray RefractionMaterial::ScatteredRay(const Ray refRay, HitRecord& rec, const double spectrum, Random& rand) const {
		double ir = 1.5;
		ir = 1.0 + (spectrum - color::MIN_WAVELENGTH) / (color::MAX_WAVELENGTH - color::MIN_WAVELENGTH) * 1.5;
		double refraction_ratio = rec.front_face ? (1.0 / ir) : ir;

		vec3 unit_direction = Normalize(refRay.dir);
		double cos_theta = fmin(Dot(-unit_direction, rec.normal), 1.0);
		double sin_theta = sqrt(1.0 - cos_theta * cos_theta);

		bool cannot_refract = refraction_ratio * sin_theta > 1.0;
		vec3 direction;

		if (cannot_refract || reflectance(cos_theta, refraction_ratio) > rand.RandomGenerate())
			direction = reflect(unit_direction, rec.normal);
		else
			direction = refract(unit_direction, rec.normal, refraction_ratio);

		auto scattered = Ray(rec.point, direction);
		return scattered;
	}

	DiffuseMaterial::DiffuseMaterial(vec3 color, vec3 emission)
		: Material(color, emission)
	{}
	Ray DiffuseMaterial::ScatteredRay(const Ray refRay, HitRecord& rec, const double spectrum, Random& rand) const {
		Ray nextRay;
		const auto orienting_normal = Dot(rec.normal, refRay.dir) < 0.0
			? rec.normal : (rec.normal * -1.0);
		nextRay.org = rec.point;
		vec3 u, v, w;
		w = orienting_normal;
		const auto r1 = 2 * PI * rand.RandomGenerate();
		const auto r2 = rand.RandomGenerate();
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
	Ray SphereObject::ScatteredRay(const Ray& refRay, HitRecord& hitrecord, const double spectrum, Random& rand) const {
		const auto distance = HitDistance(refRay);
		const auto hitPoint = refRay.dir * distance + refRay.org;
        const auto normalVec = Normalize(hitPoint - this->position());

		hitrecord.t = distance;
		hitrecord.point = refRay.org + distance * refRay.dir;
		vec3 outwardNormal = (hitrecord.point - this->position()) / this->size();
		hitrecord.set_face_normal(refRay, outwardNormal);
		//hitrecord.mat_ptr = material;

		return material->ScatteredRay(refRay, hitrecord ,spectrum, rand);
	}
	//vec3 SphereObject::Normal(vec3) const { return vec3(); }

	void SphereObject::SetPosition(vec3 posi) { m_position = posi; }
	void SphereObject::SetSize(float size) { m_size = size; }

	vec3 SphereObject::position() const { return m_position; }
	float SphereObject::size() const { return m_size; }

}