#include "SceneObject.h"
#include "Constant.h"
#include "Color.h"


namespace png {
	double kd(const vec3& color) {
		return std::max(std::max(color.x, color.y), color.z);
	}

	Material::Material(Texture* color, Texture* emission)
		: m_color(color)
		, m_emission(emission)
	{}

	vec3 Material::color(double x, double y) const { return m_color->GetColor(x, y); }
	vec3 Material::emission(double x, double y) const { return m_emission->GetColor(x, y); }

	RefractionMaterial::RefractionMaterial(const TransparentMaterialType transparentMaterialType, Texture* color, Texture* emission)
		: Material(color, emission)
		, m_transparentMaterialType(transparentMaterialType)
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
		double ir = 1.3;
		if (spectrum != -1) {
			ir = refractiveIndex(m_transparentMaterialType, spectrum);
		}
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

	DiffuseMaterial::DiffuseMaterial(Texture* color, Texture* emission)
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

		const double minValue = 1e-5;
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

		return material->ScatteredRay(refRay, hitrecord, spectrum, rand);
	}
	vec3 SphereObject::color(const vec3& point) const {
		return material->color(0,0);
	}
	vec3 SphereObject::emission(const vec3& point) const {
		return material->emission(0, 0);
	}

	void SphereObject::SetPosition(vec3 posi) { m_position = posi; }
	void SphereObject::SetSize(float size) { m_size = size; }

	vec3 SphereObject::position() const { return m_position; }
	float SphereObject::size() const { return m_size; }


	MeshObject::MeshObject(std::string fileName, Material* mat)
		: SceneObject(mat) {

	}
	MeshObject::MeshObject(std::vector<vec3> mesh, Material* mat)
		: SceneObject(mat) {
		this->mesh = mesh;
	}

	void HitDistanceAll(const std::vector<vec3>& mesh, const Ray& ray, double& distance, int& hitMeshIndex) {
		auto t = std::numeric_limits<double>::max();
		for (int meshIndex = 0; meshIndex < mesh.size() / 3; ++meshIndex) {
			auto& p0 = mesh[meshIndex * 3 + 0];
			auto& p1 = mesh[meshIndex * 3 + 1];
			auto& p2 = mesh[meshIndex * 3 + 2];

			auto triangleNormal = Cross(p0 - p2, p1 - p2) / Magnitude(Cross(p0 - p2, p1 - p2));
			auto tmp_t = Dot(triangleNormal, p2 - ray.org) / (Dot(triangleNormal, ray.dir));

			if (t > tmp_t) {

				// exist point at triangle
				auto point = ray.org + tmp_t * ray.dir;
				auto cross01 = Cross(p1 - p0, point - p0);
				auto cross12 = Cross(p2 - p1, point - p1);
				auto cross20 = Cross(p0 - p2, point - p2);

				if (Dot(cross01, cross12) > 0 && Dot(cross12, cross20) > 0) {
					t = tmp_t;
					hitMeshIndex = meshIndex;
				}
			}
		}
		if (t == std::numeric_limits<double>::max()) {
			distance = -1;
			hitMeshIndex = -1;
			return;
		}
		else {
			distance = t;
			return;
		}
	}

	double MeshObject::HitDistance(const Ray& ray) const {
		double t;
		int meshIndex;
		HitDistanceAll(mesh, ray, t, meshIndex);
		return t;
	}

	Ray MeshObject::ScatteredRay(const Ray& refRay, HitRecord& hitrecord, const double spectrum, Random& rand) const {
		double distance;
		int meshIndex;
		HitDistanceAll(mesh, refRay, distance, meshIndex);
		const auto p0 = mesh[meshIndex * 3 + 0];
		const auto p1 = mesh[meshIndex * 3 + 1];
		const auto p2 = mesh[meshIndex * 3 + 2];
		auto triangleNormal = Cross(p0 - p2, p1 - p2) / Magnitude(Cross(p0 - p2, p1 - p2));

		hitrecord.t = distance;
		hitrecord.point = refRay.org + distance * refRay.dir;
		vec3 outwardNormal = triangleNormal;
		hitrecord.set_face_normal(refRay, outwardNormal);
		//hitrecord.mat_ptr = material;

		return material->ScatteredRay(refRay, hitrecord, spectrum, rand);
	}

	vec3 MeshObject::color(const vec3& point) const {
		return material->color(0,0);
	}
	vec3 MeshObject::emission(const vec3& point) const {
		return material->emission(0,0);
	}

	BoxObject::BoxObject(const vec3& offset, Material* mat)
		: MeshObject(std::vector<vec3>({
				//front
				vec3(1,1,-1), vec3(-1,1,-1), vec3(-1,-1,-1),
				vec3(1,1,-1), vec3(-1,-1,-1), vec3(1,-1,-1),
				//top
				vec3(1,1,1), vec3(-1,1,1), vec3(-1,1,-1),
				vec3(1,1,1), vec3(-1,1,-1), vec3(1,1,-1),
				//bottum
				vec3(1,-1,-1), vec3(-1,-1,-1), vec3(-1,-1,1),
				vec3(1,-1,-1), vec3(-1,-1,1), vec3(1,-1,1),
				//right
				vec3(1,1,1), vec3(1,1,-1), vec3(1,-1,-1),
				vec3(1,1,1), vec3(1,-1,-1), vec3(1,-1,1),
				//left
				vec3(-1,1,-1), vec3(-1,1,1), vec3(-1,-1,1),
				vec3(-1,1,-1), vec3(-1,-1,1), vec3(-1,-1,-1),
				//back
				vec3(-1,1,1), vec3(1,1,1), vec3(1,-1,1),
				vec3(-1,1,1), vec3(1,-1,1), vec3(-1,-1,1),
			}), mat) {

		//offset
		for (int i = 0; i < mesh.size(); ++i) {
			mesh[i] += offset;
		}
	}

	/*
	BoxObject::BoxObject(Material* mat)
		: MeshObject(std::vector<vec3>({
		//front
		vec3(1,1,-1), vec3(-1,1,-1), vec3(-1,-1,-1),
		//vec3(1,1,-1), vec3(-1,-1,-1), vec3(1,-1,-1)
			}), mat) {

	}
	*/

	vec3 BoxObject::color(const vec3& point) const {
		return material->color(0,0);
	}
	vec3 BoxObject::emission(const vec3& point) const {
		return material->emission(0,0);
	}

	PlaneObject::PlaneObject(const vec3& offset, const vec3& up, const vec3& right, Material* mat)
	: MeshObject(std::vector<vec3>(), mat)
	, m_up(up)
	, m_right(right)
	, m_offset(offset){
		mesh.push_back(up + right + offset);
		mesh.push_back(up - right + offset);
		mesh.push_back(- up - right + offset);
		//
		mesh.push_back(up + right + offset);
		mesh.push_back(-up - right + offset);
		mesh.push_back(-up + right + offset);
	}

	//color
	vec3 PlaneObject::color(const vec3& point) const {
		auto pointInTexture = point - (-m_up - m_right);
		auto dotX = Dot(pointInTexture, m_right);
		auto dotY = Dot(pointInTexture, m_up);

		return vec3(dotX, dotY, 0);
	}
	vec3 PlaneObject::emission(const vec3& point) const {
		auto pointInTexture = point - (-m_up - m_right);
		auto dotX = Dot(pointInTexture, 2.0*m_right) / Dot(2.0*m_right, 2.0*m_right);
		auto dotY = Dot(pointInTexture, 2.0 * m_up) / Dot(2.0 * m_up, 2.0 * m_up);

		return material->emission(dotX, dotY);
	}

}