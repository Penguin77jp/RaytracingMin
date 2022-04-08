#include "SceneObject.h"
#include "Constant.h"
#include "Color.h"
#include "SettingData.h"

namespace png {
	double kd(const vec3& color) {
		return std::max(std::max(color.x, color.y), color.z);
	}
	void RayCasting(const Ray& ray, const SettingData& data, int& objectIndex, double& distance) {
		objectIndex = -1;
		distance = std::numeric_limits<double>::max();
		for (int i = 0; i < data.object.size(); ++i) {
			auto& obj = data.object[i];
			auto tmp_dis = obj->HitDistance(ray);
			if (tmp_dis < distance && tmp_dis > 0) {
				distance = tmp_dis;
				objectIndex = i;
			}
		}
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

	inline vec3 refract(const vec3& incomingDir, const vec3& normal, double n1, double n2) {
		auto cos2 = 1.0 - pow2(n1 / n2) * (1.0 - pow2(Dot(-incomingDir, normal)));
		if (cos2 <= 0.0) {
			// full refract
			//return std::nullopt;
		}
		return n1 / n2 * incomingDir - (n1 / n2 * Dot(incomingDir, normal) + std::sqrt(cos2)) * normal;
	}

	Ray RefractionMaterial::ScatteredRay(const Ray refRay, HitRecord& rec, const png::SettingData& data, double& out_pdfWeight, const double spectrum, Random& rand) const {
		out_pdfWeight = 1.0;
		double ir = refractiveIndex(m_transparentMaterialType, spectrum);
		double n1, n2;
		if (rec.front_face) {
			n1 = 1.0;
			n2 = ir;
		}
		else {
			n1 = ir;
			n2 = 1.0;
		}
		double refraction_ratio = n1 / n2;

		vec3 unit_direction = Normalize(refRay.dir);
		double cos_theta = fmin(Dot(-unit_direction, rec.normal), 1.0);
		double sin_theta = sqrt(1.0 - cos_theta * cos_theta);

		bool cannot_refract = refraction_ratio * sin_theta > 1.0;
		vec3 direction;

		if (cannot_refract || reflectance(cos_theta, refraction_ratio) > rand.RandomGenerate())
			direction = reflect(unit_direction, rec.normal);
		else
			direction = refract(unit_direction, rec.normal, n1, n2);

		auto scattered = Ray(rec.point, direction);
		return scattered;
	}

	DiffuseMaterial::DiffuseMaterial(Texture* color, Texture* emission, PDF_TYPE pdf_type)
		: Material(color, emission)
		, m_PDF_TYPE(pdf_type)
	{}
	Ray DiffuseMaterial::ScatteredRay(const Ray refRay, HitRecord& rec, const png::SettingData& data, double& out_pdfWeight, const double spectrum, Random& rand) const {
		if (m_PDF_TYPE == PDF_TYPE::CosImportSample) {
			out_pdfWeight = 1.0;
			Ray nextRay;
			const auto orienting_normal = Dot(rec.normal, refRay.dir) < 0.0
				? rec.normal : (rec.normal * -1.0);
			nextRay.org = rec.point;
			vec3 u, v, w;
			w = orienting_normal;
			if (fabs(w.x) > EPS) {
				u = Normalize(Cross(vec3(0, 1, 0), w));
			}
			else {
				u = Normalize(Cross(vec3(1, 0, 0), w));
			}
			v = Cross(w, u);
			auto xi_theta = 0.5 * std::acos(-2.0 * rand.RandomGenerate() + 1);
			auto xi_phi = 2.0 * PI * rand.RandomGenerate();
			nextRay.dir = Normalize(
				u * std::sin(xi_theta) * std::cos(xi_phi) +
				v * std::sin(xi_theta) * std::sin(xi_phi) +
				w * std::cos(xi_theta)
			);
			return nextRay;
		}
		else if (m_PDF_TYPE == PDF_TYPE::LightSample) {
			int sampleObjectIndex = -1;
			vec3 randomSurfacePoint;
			vec3 dir;
			vec3 lightNormal;
			while (true) {
				sampleObjectIndex = (int)(rand.RandomGenerate() * data.object.size());
				// sample point in selected surface
				auto& sampleObject = data.object[sampleObjectIndex];
				randomSurfacePoint = sampleObject->RandomSurfacePoint(rand);
				dir = Normalize(randomSurfacePoint - rec.point);
				lightNormal = sampleObject->Normal(randomSurfacePoint);
				if (Dot(dir, lightNormal) <= 0) {
					break;
				}
			}
			// geometory calculate
			// possible path throw from point x to surface point
			auto pointx = rec.point;
			auto distance = Magnitude(randomSurfacePoint - pointx);
			int hitObjectIndex = -1;
			double rayCastingDistance = 0.0;
			RayCasting(Ray(pointx, dir), data, hitObjectIndex, rayCastingDistance);
			if (hitObjectIndex == sampleObjectIndex && std::abs(rayCastingDistance - distance) < EPS) {
				// return pdf weight
				out_pdfWeight = std::abs(Dot(rec.normal, dir) * Dot(lightNormal, dir)) / pow2(distance);
				// return scattered ray
				return Ray(pointx, dir);
			}
			else {
				out_pdfWeight = 0.0;
				return Ray();
			}
		}
		else if (m_PDF_TYPE == PDF_TYPE::Mixture) {
			out_pdfWeight = 0.0;
			if (rand.RandomGenerate() < 0.5) {
				int sampleObjectIndex = -1;
				vec3 randomSurfacePoint;
				vec3 dir;
				vec3 lightNormal;
				while (true) {
					sampleObjectIndex = (int)(rand.RandomGenerate() * data.object.size());
					// sample point in selected surface
					auto& sampleObject = data.object[sampleObjectIndex];
					randomSurfacePoint = sampleObject->RandomSurfacePoint(rand);
					dir = Normalize(randomSurfacePoint - rec.point);
					lightNormal = sampleObject->Normal(randomSurfacePoint);
					if (Dot(dir, lightNormal) <= 0) {
						break;
					}
				}
				// geometory calculate
				// possible path throw from point x to surface point
				auto pointx = rec.point;
				auto distance = Magnitude(randomSurfacePoint - pointx);
				int hitObjectIndex = -1;
				double rayCastingDistance = 0.0;
				RayCasting(Ray(pointx, dir), data, hitObjectIndex, rayCastingDistance);
				if (hitObjectIndex == sampleObjectIndex && std::abs(rayCastingDistance - distance) < EPS) {
					// return pdf weight
					out_pdfWeight = std::abs(Dot(rec.normal, dir) * Dot(lightNormal, dir)) / pow2(distance);
					// return scattered ray
					return Ray(pointx, dir);
				}
				else {
					out_pdfWeight = 0.0;
					return Ray();
				}
			}
			else {
				out_pdfWeight = 1.0;
				Ray nextRay;
				const auto orienting_normal = Dot(rec.normal, refRay.dir) < 0.0
					? rec.normal : (rec.normal * -1.0);
				nextRay.org = rec.point;
				vec3 u, v, w;
				w = orienting_normal;
				if (fabs(w.x) > std::numeric_limits<float>::min()) {
					u = Normalize(Cross(vec3(0, 1, 0), w));
				}
				else {
					u = Normalize(Cross(vec3(1, 0, 0), w));
				}
				v = Cross(w, u);
				auto xi_theta = 0.5 * std::acos(-2.0 * rand.RandomGenerate() + 1);
				auto xi_phi = 2.0 * PI * rand.RandomGenerate();
				nextRay.dir = Normalize(
					u * std::sin(xi_theta) * std::cos(xi_phi) +
					v * std::sin(xi_theta) * std::sin(xi_phi) +
					w * std::cos(xi_theta)
				);
				return nextRay;
			}
		}

		return Ray();
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

		if (t1 < EPS && t2 < EPS)
			return 0;

		if (t1 > EPS) {
			return t1;
		}
		else {
			return t2;
		}
	}
	Ray SphereObject::ScatteredRay(const Ray& refRay, HitRecord& hitrecord, const SettingData& data, double& out_pdfWeight, const double spectrum, Random& rand) const {
		vec3 outwardNormal = (hitrecord.point - this->position()) / this->size();
		hitrecord.set_face_normal(refRay, outwardNormal);

		return material->ScatteredRay(refRay, hitrecord, data, out_pdfWeight, spectrum, rand);
	}
	vec3 SphereObject::RandomSurfacePoint(Random& rand) const {
		double xi_theta = rand.RandomGenerate() * PI;
		double xi_phi = 2.0 * rand.RandomGenerate() * PI;
		auto dir = vec3(sin(xi_theta) * cos(xi_phi), sin(xi_theta) * sin(xi_phi), cos(xi_theta));
		return size() * dir + position();
	}
	vec3 SphereObject::Normal(const vec3& point) const {
		return Normalize(point - position());
	}
	vec3 SphereObject::color(const vec3& point) const {
		return material->color(0, 0);
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
		this->m_mesh = mesh;

		// calculate area
		auto sumArea = 0.0;
		for (int i = 0; i < m_mesh.size() / 3; ++i) {
			auto cal = Magnitude(Cross(m_mesh[i * 3 + 1] - m_mesh[i * 3], m_mesh[i * 3 + 2] - m_mesh[i * 3]));
			m_area.push_back(cal);
			sumArea += cal;
		}
		m_SumArea = sumArea;
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

				if (Dot(cross01, cross12) >= 0 && Dot(cross12, cross20) >= 0) {
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
		HitDistanceAll(m_mesh, ray, t, meshIndex);
		return t;
	}

	Ray MeshObject::ScatteredRay(const Ray& refRay, HitRecord& hitrecord, const SettingData& data, double& out_pdfWeight, const double spectrum, Random& rand) const {
		double distance;
		int meshIndex;
		HitDistanceAll(m_mesh, refRay, distance, meshIndex);
		const auto p0 = m_mesh[meshIndex * 3 + 0];
		const auto p1 = m_mesh[meshIndex * 3 + 1];
		const auto p2 = m_mesh[meshIndex * 3 + 2];
		auto triangleNormal = Cross(p0 - p2, p1 - p2) / Magnitude(Cross(p0 - p2, p1 - p2));

		vec3 outwardNormal = triangleNormal;
		hitrecord.set_face_normal(refRay, outwardNormal);
		//hitrecord.mat_ptr = material;

		return material->ScatteredRay(refRay, hitrecord, data, out_pdfWeight, spectrum, rand);
	}
	double SumArea(const std::vector<double> elements, const int index) {
		double sum = 0.0;
		for (int i = 0; i <= index; ++i) {
			sum += elements[i];
		}
		return sum;
	}
	vec3 Square2Triangle(vec3 point) {
		vec3 out;
		if (point.y > point.x) {
			out.x = point.x * 0.5;
			out.y = point.y - out.x;
		}
		else {
			out.y = point.y * 0.5;
			out.x = point.x - out.y;
		}
		return out;
	}
	vec3 MeshObject::RandomSurfacePoint(Random& rand) const {
		int index = -1;
		{
			double randomIndex = rand.RandomGenerate();
			for (int i = 0; i < m_area.size(); ++i) {
				if (randomIndex < SumArea(m_area, i) / m_SumArea) {
					index = i;
					break;
				}
			}
		}
		vec3 randomPoint;
		{
			vec3 randomSquare = vec3(rand.RandomGenerate(), rand.RandomGenerate(), 0);
			vec3 randomTriangle = Square2Triangle(randomSquare);
			randomPoint = (m_mesh[index * 3 + 1] - m_mesh[index * 3]) * randomTriangle.x
				+ (m_mesh[index * 3 + 2] - m_mesh[index * 3]) * randomTriangle.y
				+ m_mesh[index * 3];
		}

		return randomPoint;
	}
	vec3 MeshObject::Normal(const vec3& point) const {
		for (int i = 0; i < m_mesh.size() / 3; ++i) {
			auto p0 = m_mesh[i * 3] - point;
			auto p1 = m_mesh[i * 3 + 1] - point;
			auto p2 = m_mesh[i * 3 + 2] - point;
			auto cross01 = Cross(p0, p1);
			auto cross12 = Cross(p1, p2);
			auto cross20 = Cross(p2, p0);
			if (Dot(cross01, cross12) >= 0 && Dot(cross12, cross20) >= 0) {
				return Normalize(Cross(m_mesh[i*3+1]-m_mesh[i*3], m_mesh[i*3+2]-m_mesh[i*3]));
			}
		}

		return vec3();
	}

	vec3 MeshObject::color(const vec3& point) const {
		return material->color(0, 0);
	}
	vec3 MeshObject::emission(const vec3& point) const {
		return material->emission(0, 0);
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
		for (int i = 0; i < m_mesh.size(); ++i) {
			m_mesh[i] += offset;
		}
	}
	vec3 BoxObject::RandomSurfacePoint(Random& rand) const {
		return vec3();
	}
	vec3 BoxObject::Normal(const vec3& point) const {
		return vec3();
	}

	vec3 BoxObject::color(const vec3& point) const {
		return material->color(0, 0);
	}
	vec3 BoxObject::emission(const vec3& point) const {
		return material->emission(0, 0);
	}

	PlaneObject::PlaneObject(const vec3& offset, const vec3& up, const vec3& right, Material* mat)
		: MeshObject(std::vector<vec3>(), mat)
		, m_up(up)
		, m_right(right)
		, m_offset(offset) {
		m_mesh.push_back(up + right + offset);
		m_mesh.push_back(up - right + offset);
		m_mesh.push_back(-up - right + offset);
		//
		m_mesh.push_back(up + right + offset);
		m_mesh.push_back(-up - right + offset);
		m_mesh.push_back(-up + right + offset);
	}
	vec3 PlaneObject::RandomSurfacePoint(Random& rand) const {
		return vec3();
	}
	vec3 PlaneObject::Normal(const vec3& point) const {
		return vec3();
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
		auto dotX = Dot(pointInTexture, 2.0 * m_right) / Dot(2.0 * m_right, 2.0 * m_right);
		auto dotY = Dot(pointInTexture, 2.0 * m_up) / Dot(2.0 * m_up, 2.0 * m_up);

		return material->emission(dotX, dotY);
	}

}