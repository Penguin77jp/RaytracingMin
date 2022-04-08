#include "Camera.h"
#include "Ray.h"
#include "Random.h"
#include "color.h"
#include "Transparent.h"

#include <numbers>

// debug
#include <fstream>

bool inline DEBUG() {
#ifdef _DEBUG
	return true;
#else
	return false;
#endif
}

namespace png {
	void kai(double a, double b, double c, std::optional<double>& t1, std::optional<double>& t2) {
		double b24ac = b * b - 4 * a * c;
		double eps = 1e-10;
		if (eps < b24ac) {
			t1 = (-b + sqrt(b24ac)) / (2 * a);
			t2 = (-b - sqrt(b24ac)) / (2 * a);
			t1 = t1 > 0.0 ? t1 : std::nullopt;
			t2 = t2 > 0.0 ? t2 : std::nullopt;
		}
		else if (-eps < b24ac && b24ac < eps) {
			t1 = (-b) / (2 * a);
			t1 = t1 > 0.0 ? t1 : std::nullopt;
			t2 = std::nullopt;
		}
		else {
			t1 = std::nullopt;
			t2 = std::nullopt;
		}
	}
	std::optional<double> hitSphere(const vec3& spherePosition, const double sphereRadius, const Ray& ray) {
		std::optional<double> t1, t2;
		kai(Dot(ray.dir, ray.dir), -2.0 * Dot(ray.dir, spherePosition - ray.org), Dot(spherePosition - ray.org, spherePosition - ray.org) - sphereRadius * sphereRadius, t1, t2);
		if (t1.has_value() && t2.has_value()) {
			return t2;
		}
		else if (t1.has_value()) {
			return t1;
		}
		else {
			return std::nullopt;
		}
	}

	std::optional<png::Ray> Lens(const LensData& lensData, const vec3& aperturePosi, const double apertureRadius, png::Ray ray, const double spectrum) {
		const auto onlyPathedRay = true;
		std::ofstream file;
		if (DEBUG())
			file = std::ofstream("ray.txt", std::ios_base::app);
		auto rayCal = ray;
		// first surface
		{
			auto sphperePosition = lensData.posi + vec3(0, 0, lensData.radius_front - lensData.thickness * 0.5);
			auto sphereRadius = std::abs(lensData.radius_front);
			auto sphereDistance = hitSphere(sphperePosition, sphereRadius, rayCal);
			if (!sphereDistance.has_value()) {
				if (DEBUG() && !onlyPathedRay) {
					file << std::string(rayCal.org) << std::endl;
					file << std::string(rayCal.org + rayCal.dir) << std::endl;
				}
				return std::nullopt;
			}
			// check aperture
			{
				double apertureDistance = (aperturePosi.z - rayCal.org.z) / rayCal.dir.z;
				if (apertureDistance <= sphereDistance.value()) {
					auto point = ray.org + apertureDistance * ray.dir;
					auto r_pow2 = point.x * point.x + point.y * point.y;
					if (r_pow2 > apertureRadius * apertureRadius) {
						if (DEBUG() && !onlyPathedRay) {
							file << std::string(rayCal.org) << std::endl;
							file << std::string(point) << std::endl;
						}
						return std::nullopt;
					}
				}
			}
			// refractive
			auto eta_over_eta = 1.0 / refractiveIndex(lensData.traType, spectrum);
			auto sphereHitpoint = rayCal.org + sphereDistance.value() * rayCal.dir;
			auto normalVec = Normalize(sphereHitpoint - sphperePosition);
			normalVec = Dot(rayCal.dir, normalVec) > 0 ? normalVec : -normalVec;
			if (DEBUG())
				file << std::string(rayCal.org) << std::endl;
			rayCal.org = sphereHitpoint;
			auto outDir_root = std::sqrt(1.0 - pow2(eta_over_eta) * (1.0 - pow2(Dot(normalVec, rayCal.dir))));
			rayCal.dir = outDir_root * normalVec + eta_over_eta * (rayCal.dir - Dot(normalVec, rayCal.dir) * normalVec);
			if (DEBUG())
				file << std::string(rayCal.org) << std::endl;
			/*
			auto sphereHitpoint = rayCal.org + sphereDistance.value() * rayCal.dir;
			vec3 outcomingVec;
			{
				auto normalVec = Normalize(sphereHitpoint - sphperePosition);
				vec3 w, u, v;
				w = Dot(normalVec, rayCal.dir) < 0.0 ? normalVec : (-1.0 * normalVec);
				if (abs(w.x) > 1e-10) {
					u = Normalize(Cross(vec3(0.0, 1.0, 0.0), w));
				}
				else {
					u = Normalize(Cross(vec3(1.0, 0.0, 0.0), w));
				}
				v = Cross(w, u);
				vec3 outcomingVecU, outcomingVecV;
				{
					vec3 incomingVec = u * Dot(-1.0 * rayCal.dir, u);
					double incomingVec_length = Magnitude(incomingVec);
					vec3 incomingVec_normal = incomingVec / incomingVec_length;
					auto T_root = 1.0 - pow2(refractiveIndex(lensData.traType, spectrum)) * (1.0 - pow2(Dot(incomingVec_normal, w)));
					// zen hansya
					if (T_root < 0) {
						return std::nullopt;
					}
					outcomingVecU = -std::sqrt(T_root) * w - refractiveIndex(lensData.traType, spectrum) * (incomingVec_normal - w * Dot(incomingVec_normal, w));

				}
				{
					vec3 incomingVec = v * Dot(-1.0 * rayCal.dir, v);
					double incomingVec_length = Magnitude(incomingVec);
					vec3 incomingVec_normal = incomingVec / incomingVec_length;
					auto T_root = 1.0 - pow2(refractiveIndex(lensData.traType, spectrum)) * (1.0 - pow2(Dot(incomingVec_normal, w)));
					// zen hansya
					if (T_root < 0) {
						return std::nullopt;
					}
					outcomingVecV = -std::sqrt(T_root) * w - refractiveIndex(lensData.traType, spectrum) * (incomingVec_normal - w * Dot(incomingVec_normal, w));
					outcomingVecV = outcomingVecV * incomingVec_length;
				}
				auto hoge = outcomingVecU + outcomingVecV;
				auto fuga = Magnitude(hoge);
				outcomingVec = Normalize(outcomingVecU + outcomingVecV);
			}
			rayCal = Ray(sphereHitpoint, outcomingVec);
			*/
		}

		// second surface
		{
			auto sphperePosition = lensData.posi - vec3(0, 0, lensData.radius_rear - lensData.thickness * 0.5);
			auto sphereRadius = std::abs(lensData.radius_rear);
			auto sphereDistance = hitSphere(sphperePosition, sphereRadius, rayCal);
			if (!sphereDistance.has_value()) {
				if (DEBUG() && !onlyPathedRay) {
					file << std::string(rayCal.org) << std::endl;
					file << std::string(rayCal.org + rayCal.dir) << std::endl;
				}
				return std::nullopt;
			}
			// check aperture
			{
				double apertureDistance = (aperturePosi.z - rayCal.org.z) / rayCal.dir.z;
				if (apertureDistance <= sphereDistance.value()) {
					auto point = ray.org + apertureDistance * ray.dir;
					auto r_pow2 = point.x * point.x + point.y * point.y;
					if (r_pow2 > apertureRadius * apertureRadius) {
						if (DEBUG() && !onlyPathedRay) {
							file << std::string(rayCal.org) << std::endl;
							file << std::string(point) << std::endl;
						}
						return std::nullopt;
					}
				}
			}
			// refractive
			auto eta_over_eta = refractiveIndex(lensData.traType, spectrum);
			auto sphereHitpoint = rayCal.org + sphereDistance.value() * rayCal.dir;
			auto normalVec = Normalize(sphereHitpoint - sphperePosition);
			normalVec = Dot(rayCal.dir, normalVec) > 0 ? normalVec : -normalVec;
			if (DEBUG())
				file << std::string(rayCal.org) << std::endl;
			rayCal.org = sphereHitpoint;
			auto outDir_root = std::sqrt(1.0 - pow2(eta_over_eta) * (1.0 - pow2(Dot(normalVec, rayCal.dir))));
			rayCal.dir = outDir_root * normalVec + eta_over_eta * (rayCal.dir - Dot(normalVec, rayCal.dir) * normalVec);
			if (DEBUG()) {
				file << std::string(rayCal.org) << std::endl;
				file << std::string(rayCal.org) << std::endl;
				file << std::string(rayCal.org + rayCal.dir) << std::endl;
			}


			/*
			auto sphereHitpoint = rayCal.org + sphereDistance.value() * rayCal.dir;
			vec3 outcomingVec;
			{
				auto normalVec = Normalize(sphereHitpoint - sphperePosition);
				vec3 w, u, v;
				w = Dot(normalVec, rayCal.dir) < 0.0 ? normalVec : (-1.0 * normalVec);
				if (abs(w.x) > 1e-10) {
					u = Normalize(Cross(vec3(0.0, 1.0, 0.0), w));
				}
				else {
					u = Normalize(Cross(vec3(1.0, 0.0, 0.0), w));
				}
				v = Cross(w, u);
				vec3 outcomingVecU, outcomingVecV;
				{
					vec3 incomingVec = Normalize(u * Dot(-1.0 * rayCal.dir, u));
					auto T_root = 1.0 - pow2(refractiveIndex(lensData.traType, spectrum)) * (1.0 - pow2(Dot(incomingVec, w)));
					// zen hansya
					if (T_root < 0) {
						return std::nullopt;
					}
					outcomingVecU = -std::sqrt(T_root) * w - refractiveIndex(lensData.traType, spectrum) * (incomingVec - w * Dot(incomingVec, w));
				}
				{
					vec3 incomingVec = Normalize(v * Dot(-1.0 * rayCal.dir, v));
					auto T_root = 1.0 - pow2(refractiveIndex(lensData.traType, spectrum)) * (1.0 - pow2(Dot(incomingVec, w)));
					// zen hansya
					if (T_root < 0) {
						return std::nullopt;
					}
					outcomingVecV = -std::sqrt(T_root) * w - refractiveIndex(lensData.traType, spectrum) * (incomingVec - w * Dot(incomingVec, w));
				}
				outcomingVec = Normalize(outcomingVecU + outcomingVecV);
			}
			rayCal = Ray(sphereHitpoint, outcomingVec);
			*/
		}

		// check aperture
		{
			double apertureDistance = (aperturePosi.z - rayCal.org.z) / rayCal.dir.z;
			auto point = ray.org + apertureDistance * ray.dir;
			auto r_pow2 = point.x * point.x + point.y * point.y;
			if (r_pow2 > apertureRadius * apertureRadius) {
				if (DEBUG() && !onlyPathedRay) {
					file << std::string(rayCal.org) << std::endl;
					file << std::string(point) << std::endl;
				}
				return std::nullopt;
			}
		}

		return rayCal;
	}

	Camera::Camera(const SettingData& data)
		:m_data(data),
		m_target(data.cameraTarget),
		m_origin(data.cameraOrigin),
		m_fov(data.fov),
		m_superSamples(data.superSamples)
	{}

	vec3 Camera::target() const { return m_target; }
	vec3 Camera::origin() const { return m_origin; }
	double Camera::fov() const { return m_fov; }
	int Camera::superSamples() const { return m_superSamples; }

	NoLensCamera::NoLensCamera(SettingData& data)
		: Camera(data)
	{
		vec3 upVec = vec3(0, 1, 0);
		m_direction = Normalize(target() - origin());
		m_camX = -Normalize(Cross(m_direction, upVec));
		m_camY = Cross(m_camX, m_direction);
		m_camZ = m_direction;
		//fov
		m_fovx = fov();
		m_fovy = m_fovx * data.height / data.width;
	}
	void NoLensCamera::GenerateRay(int x, int y, int superSampleX, int superSampleY, double spectrum, Ray& rayIncomingSensor, std::optional<Ray>& generatedRay, Random& rand) {
		rayIncomingSensor = Ray();
		rayIncomingSensor.dir = Normalize(target() - origin());
		const auto rate = 1.0 / (1 + superSamples());
		vec3 dir = Normalize(
			m_camX * m_fovx * (2.0f * ((double)x + rate * superSampleX) / m_data.width - 1.0f) +
			m_camY * m_fovy * (2.0f * ((double)y + rate * superSampleY) / m_data.height - 1.0f) +
			m_camZ
		);
		generatedRay = Ray(m_data.cameraOrigin, dir);
	}

	LensCamera::LensCamera(const std::vector<LensData>& lensData, const SettingData& data)
		: m_lensData(lensData)
		, Camera(data)
	{
		vec3 upVec = vec3(0, 1, 0);
		m_direction = Normalize(target() - origin());
		m_camX = -Normalize(Cross(m_direction, upVec));
		m_camY = Cross(m_camX, m_direction);
		m_camZ = m_direction;
		//fov
		m_fovx = fov();
		m_fovy = m_fovx * data.height / data.width;
	}

	LensCamera::LensCamera(const SettingData& data, const PrototypeLensType prototypeType)
		: Camera(data)
	{
		vec3 upVec = vec3(0, 1, 0);
		m_direction = Normalize(target() - origin());
		m_camX = -Normalize(Cross(m_direction, upVec));
		m_camY = Cross(m_camX, m_direction);
		m_camZ = m_direction;
		//fov
		m_fovx = fov();
		m_fovy = m_fovx * data.height / data.width;

		if (prototypeType == PrototypeLensType::Concave) {
			m_aperturePosition = vec3(0, 0, 0.6);
			m_apertureRadius = 0.5;
			auto lensData1 = LensData();
			lensData1.posi = vec3(0, 0, 0.6);
			lensData1.thickness = 1.0;
			lensData1.radius_front = -10.0;
			lensData1.radius_rear = 1e10;
			lensData1.traType = TransparentMaterialType::HighVariance;
			m_lensData = std::vector<LensData>{ lensData1 };
		}
		else if (prototypeType == PrototypeLensType::Convex) {
			m_aperturePosition = vec3(0, 0, 0.6);
			m_apertureRadius = 1.0;
			auto lensData1 = LensData();
			lensData1.posi = vec3(0, 0, 0.6);
			lensData1.thickness = 1.0;
			lensData1.radius_front = 2.0;
			lensData1.radius_rear = 2.0;
			lensData1.traType = TransparentMaterialType::BK7;
			m_lensData = std::vector<LensData>{ lensData1 };
		}
	}

	void LensCamera::GenerateRay(int x, int y, int superSampleX, int superSampleY, double spectrum, Ray& rayIncomingSensor, std::optional<Ray>& generatedRay, Random& rand) {
		Ray lensRay;
		{
			const auto rate = 1.0 / (1 + superSamples());
			const auto maxImageSensorResolution = std::max(m_data.width, m_data.height);
			const auto indexX = (double)x + rate * superSampleX;
			const auto indexY = (double)y + rate * superSampleY;
			lensRay.org = vec3((2.0 * indexX / m_data.width - 1.0) * m_data.width / maxImageSensorResolution
				, (2.0 * indexY / m_data.height - 1.0) * m_data.height / maxImageSensorResolution
				, 0.0);
			std::optional<Ray> outRay;
			//double xi_theta = 0.5 * std::acos(-2.0 * rand.RandomGenerate() + 1.0);
			//double xi_phi = 2.0 * std::numbers::pi * rand.RandomGenerate();
			//lensRay.dir = Normalize(vec3(std::sin(xi_theta) * std::cos(xi_phi), std::sin(xi_theta) * std::sin(xi_phi), std::cos(xi_theta)));

			//debug
			while (true) {
				double xi_theta = 0.5 * std::acos(-2.0 * rand.RandomGenerate() + 1.0);
				double xi_phi = 2.0 * std::numbers::pi * rand.RandomGenerate();
				lensRay.dir = Normalize(vec3(std::sin(xi_theta) * std::cos(xi_phi), std::sin(xi_theta) * std::sin(xi_phi), std::cos(xi_theta)));
				if (Dot(lensRay.dir, vec3(0, 0, 1)) > 1.0 - 1e-1) {
					break;
				}
			}

			rayIncomingSensor = lensRay;
		}

		for (int i = 0; i < m_lensData.size(); ++i) {
			auto resultRay = Lens(m_lensData[i], m_aperturePosition, m_apertureRadius, lensRay, spectrum);
			if (!resultRay.has_value()) {
				generatedRay = std::nullopt;
				return;
			}
			lensRay = resultRay.value();
		}
		Ray generatedRayWorldSpace;
		generatedRayWorldSpace.org = m_camX * m_fovx * lensRay.org.x + m_camY * m_fovy * lensRay.org.y + m_data.cameraOrigin;
		generatedRayWorldSpace.dir = Normalize(m_camX * m_fovx * lensRay.dir.x + m_camY * m_fovy * lensRay.dir.y + m_camZ * lensRay.dir.z);
		generatedRay = generatedRayWorldSpace;
	}

	LensSystem::LensSystem(const std::string fileName, const SettingData& data)
		: Camera(data) {
		vec3 upVec = vec3(0, 1, 0);
		m_direction = Normalize(target() - origin());
		m_camX = -Normalize(Cross(m_direction, upVec));
		m_camY = Cross(m_camX, m_direction);
		m_camZ = m_direction;
		m_targetDistance = Magnitude(target() - origin());
		//fov
		m_fovx = fov();
		m_fovy = m_fovx * data.height / data.width;

		std::ifstream file(fileName);
		nlohmann::json json;
		file >> json;

		for (auto& it : json.items()) {
			png::LensSurface tmp_lens;
			for (auto& itLens : it.value().items()) {
				if (itLens.key() == "r") {
					tmp_lens.r = itLens.value();
				}
				else if (itLens.key() == "h") {
					tmp_lens.h = itLens.value();
				}
				else if (itLens.key() == "d") {
					tmp_lens.d = itLens.value();
				}
				else if (itLens.key() == "ior") {
					tmp_lens.ior = itLens.value();
				}
			}
			m_lensSurface.push_back(tmp_lens);
		}

		// init z
		double z = 0.0;
		for (int i = m_lensSurface.size() - 1; i >= 0; --i) {
			z -= m_lensSurface[i].d;
			m_lensSurface[i].z = z;
		}

		// principle point
		for (int i = 0; i < m_lensSurface.size(); ++i) {
			auto& targetLens = m_lensSurface[i];

			if (std::abs(targetLens.r) > 0.0) {
				auto sign = targetLens.r / std::abs(targetLens.r);
				targetLens.pp = targetLens.r - sign * std::sqrt(pow2(targetLens.r) - pow2(targetLens.h));
			}
			else {
				targetLens.pp = std::nan("");
			}
		}

		// focal length
		for (int i = 0; i < m_lensSurface.size() - 1; ++i) {
			auto& targetLens = m_lensSurface[i];
			auto& nextLens = m_lensSurface[i + 1];

			if (targetLens.ior == 1.0) {
				targetLens.f = 0;
			}
			else {
				// te : thickness of lens at edge
				//auto te = targetLens.d - std::max(targetLens.pp, 0.0) - std::max(-nextLens.pp, 0.0);
				auto FLover = (targetLens.ior - 1.0) * (1.0 / targetLens.r - 1.0 / nextLens.r
					+ (targetLens.ior - 1.0) * targetLens.d / targetLens.ior / targetLens.r / nextLens.r);
				targetLens.f = 1.0 / FLover;
			}
		}
	}


	void Intersect(png::vec3 sphereCenter, LensSurface& lensData, png::Ray ray, std::optional<png::vec3>& hitpoint, png::vec3& normal) {
		if (std::abs(lensData.r) > 0.0) {
			auto p_o = sphereCenter - ray.org;
			auto b = Dot(p_o, ray.dir);
			double d4 = b * b - Dot(p_o, p_o) + pow2(lensData.r);

			if (d4 < 0.0) {
				hitpoint = std::nullopt;
				return;
			}

			double d4sqrt = std::sqrt(d4);
			double t1 = b - d4sqrt;
			double t2 = b + d4sqrt;

			double t = 0.0;
			if (t1 > png::EPS) {
				t = t1;
			}
			else if (t2 > png::EPS) {
				t = t2;
			}
			else {
				hitpoint = std::nullopt;
				return;
			}

			hitpoint = ray.org + ray.dir * t;
			if (pow2(hitpoint.value().x) + pow2(hitpoint.value().y) > pow2(lensData.h)) {
				hitpoint = std::nullopt;
				return;
			}

			normal = Normalize(hitpoint.value() - sphereCenter);
			if (Dot(ray.dir, normal) > 0) {
				normal = -normal;
			}
		}
		else {
			double t = (lensData.z - ray.org.z) / ray.dir.z;
			hitpoint = ray.org + t * ray.dir;
			if (pow2(hitpoint.value().x) + pow2(hitpoint.value().y) > pow2(lensData.h)) {
				hitpoint = std::nullopt;
				return;
			}

			normal = png::vec3(0, 0, -1);
			if (Dot(ray.dir, normal) > 0) {
				normal = -normal;
			}
		}
	}

	std::optional<png::vec3> Refract(png::vec3 incomingDir, png::vec3 lensNormal, double n1, double n2) {
		auto cos2 = 1.0 - pow2(n1 / n2) * (1.0 - pow2(Dot(-incomingDir, lensNormal)));
		if (cos2 <= 0.0) {
			// full refract
			return std::nullopt;
		}
		return n1 / n2 * incomingDir - (n1 / n2 * Dot(incomingDir, lensNormal) + std::sqrt(cos2)) * lensNormal;
	}

	std::vector<png::vec3> Pathtracing(std::vector<LensSurface>& lensData, png::Ray ray) {
		std::vector<png::vec3> points;
		double n1 = 1.0;
		points.push_back(ray.org);
		for (int i = 0; i < lensData.size(); ++i) {
			auto& targetLens = lensData[i];

			std::optional<png::vec3> hitpoint;
			png::vec3 lensNormal;
			Intersect(png::vec3(0, 0, targetLens.z + targetLens.r), targetLens, ray, hitpoint, lensNormal);
			if (!hitpoint.has_value()) {
				return std::vector<png::vec3>();
			}
			points.push_back(hitpoint.value());

			auto refractedDir = Refract(ray.dir, lensNormal, n1, targetLens.ior);

			if (!refractedDir.has_value()) {
				return std::vector<png::vec3>();
			}
			ray.dir = refractedDir.value();
			ray.org = hitpoint.value();
			n1 = targetLens.ior;
		}
		points.push_back(ray.org + ray.dir);

		return points;
	}

	void LensSystem::GenerateRay(int x, int y, int superSampleX, int superSampleY, double spectrum, Ray& rayIncomingSensor, std::optional<Ray>& generatedRay, Random& rand) {
		Ray lensRay;
		double z = 0.0;
		{
			// assumption thin lens
			auto yover = 1.0 / m_lensSurface[0].f - 1.0 / m_targetDistance;
			z = 1.0 / yover;

			/*
			* assumption thick lens
			auto Si = -std::min(m_lensSurface[1].pp, 0.0) + m_targetDistance;
			auto So = m_lensSurface[0].f * Si / (Si - m_lensSurface[0].f);
			z = So - std::max(m_lensSurface[0].pp, 0.0);
			*/
		}
		while (true) {
			{
				const auto rate = 1.0 / (1 + superSamples());
				const auto maxImageSensorResolution = std::max(m_data.width, m_data.height);
				const auto indexX = (double)x + rate * superSampleX;
				const auto indexY = (double)y + rate * superSampleY;
				lensRay.org = vec3((2.0 * indexX / m_data.width - 1.0) * m_data.width / maxImageSensorResolution
					, (2.0 * indexY / m_data.height - 1.0) * m_data.height / maxImageSensorResolution
					, m_lensSurface[0].z - z);
				std::optional<Ray> outRay;
				double xi_theta = std::acos(std::sqrt(rand.RandomGenerate()));
				double xi_phi = 2.0 * std::numbers::pi * rand.RandomGenerate();
				lensRay.dir = Normalize(vec3(std::sin(xi_theta) * std::cos(xi_phi), std::sin(xi_theta) * std::sin(xi_phi), std::cos(xi_theta)));

				rayIncomingSensor = lensRay;
			}

			{
				double n1 = 1.0;
				for (int i = 0; i < m_lensSurface.size(); ++i) {
					auto& targetLens = m_lensSurface[i];

					std::optional<png::vec3> hitpoint;
					png::vec3 lensNormal;
					Intersect(png::vec3(0, 0, targetLens.z + targetLens.r), targetLens, lensRay, hitpoint, lensNormal);
					if (!hitpoint.has_value()) {
						// no ray
						generatedRay = std::nullopt;
						continue;
						return;
					}

					auto refractedDir = Refract(lensRay.dir, lensNormal, n1, targetLens.ior);

					if (!refractedDir.has_value()) {
						// no ray
						generatedRay = std::nullopt;
						continue;
						return;
					}
					lensRay.dir = refractedDir.value();
					lensRay.org = hitpoint.value();
					n1 = targetLens.ior;
				}
			}
			break;
		}

		auto raySize = 1.0;
		Ray generatedRayWorldSpace;
		generatedRayWorldSpace.org = m_camX * raySize * lensRay.org.x + m_camY * raySize * lensRay.org.y + m_data.cameraOrigin;
		generatedRayWorldSpace.dir = Normalize(m_camX * raySize * lensRay.dir.x + m_camY * raySize * lensRay.dir.y + m_camZ * lensRay.dir.z);
		generatedRay = generatedRayWorldSpace;
	}

	void Output(std::vector<vec3>& points, std::ofstream& output) {
		for (int i = 0; i < points.size() - 1; ++i) {
			output << points[i].x << ", " << points[i].y << ", " << points[i].z << std::endl;
			output << points[i + 1].x << ", " << points[i + 1].y << ", " << points[i + 1].z << std::endl;
		}
	}

	void LensSystem::Test(int x, int y, double z, double spectrum, bool twoD, Ray& rayIncomingSensor, std::optional<Ray>& generatedRay, Random& rand, std::ofstream& output) {
		if (twoD) {

			// approxily
			auto yover = 1.0 / m_lensSurface[0].f - 1.0 / m_targetDistance;
			z = 1.0 / yover;

			// exactlly
			auto Si = -std::min(m_lensSurface[1].pp, 0.0) + m_targetDistance;
			auto So = m_lensSurface[0].f * Si / (Si - m_lensSurface[0].f);
			//z = So - std::max(m_lensSurface[0].pp, 0.0);


			std::vector<png::vec3> points;
			Ray lensRay;
			{
				lensRay.org = vec3(1.0
					//lensRay.org = vec3(2.0 * rand.RandomGenerate()-1.0
					, 0.0
					, m_lensSurface[0].z - z);
				std::optional<Ray> outRay;
				lensRay.dir = Normalize(vec3(2.0 * rand.RandomGenerate()-1.0,
					0,
					1
				));

				rayIncomingSensor = lensRay;
			}

			points.push_back(lensRay.org);

			{
				double n1 = 1.0;
				for (int i = 0; i < m_lensSurface.size(); ++i) {
					auto& targetLens = m_lensSurface[i];

					std::optional<png::vec3> hitpoint;
					png::vec3 lensNormal;
					Intersect(png::vec3(0, 0, targetLens.z + targetLens.r), targetLens, lensRay, hitpoint, lensNormal);
					if (!hitpoint.has_value()) {
						// no ray
						generatedRay = std::nullopt;
						//points.push_back(lensRay.org + lensRay.dir);
						//Output(points, output);
						points.clear();
						return;
					}

					auto refractedDir = Refract(lensRay.dir, lensNormal, n1, targetLens.ior);

					if (!refractedDir.has_value()) {
						// no ray
						generatedRay = std::nullopt;
						//Output(points, output);
						points.clear();
						return;
					}
					lensRay.dir = refractedDir.value();
					lensRay.org = hitpoint.value();
					points.push_back(hitpoint.value());
					n1 = targetLens.ior;
				}
			}
			//auto t = (0.0 - lensRay.org.x) / lensRay.dir.x;
			auto t = 100.0;
			points.push_back(lensRay.org + std::max(t, 1.0) * lensRay.dir);
			Output(points, output);


			//Output(points, output);
			//return;

			m_fovx = m_fovy = 1.0;
			Ray generatedRayWorldSpace;
			generatedRayWorldSpace.org = m_camX * m_fovx * lensRay.org.x + m_camY * m_fovy * lensRay.org.y + m_data.cameraOrigin;
			generatedRayWorldSpace.dir = Normalize(m_camX * m_fovx * lensRay.dir.x + m_camY * m_fovy * lensRay.dir.y + m_camZ * lensRay.dir.z);
			std::vector<vec3> finalPoints;
			finalPoints.push_back(generatedRayWorldSpace.org);
			finalPoints.push_back(generatedRayWorldSpace.org + 100.0 * generatedRayWorldSpace.dir);
			//Output(finalPoints, output);
			generatedRay = generatedRayWorldSpace;
		}
		else {
			std::vector<png::vec3> points;
			Ray lensRay;
			{
				const auto maxImageSensorResolution = std::max(m_data.width, m_data.height);
				const auto indexX = (double)x;
				const auto indexY = (double)y;
				lensRay.org = vec3((2.0 * indexX / m_data.width - 1.0) * m_data.width / maxImageSensorResolution
					, (2.0 * indexY / m_data.height - 1.0) * m_data.height / maxImageSensorResolution
					, m_lensSurface[0].z + z);
				std::optional<Ray> outRay;
				double xi_theta = std::acos(std::sqrt(rand.RandomGenerate()));
				double xi_phi = 2.0 * std::numbers::pi * rand.RandomGenerate();
				lensRay.dir = Normalize(vec3(std::sin(xi_theta) * std::cos(xi_phi), std::sin(xi_theta) * std::sin(xi_phi), std::cos(xi_theta)));

				rayIncomingSensor = lensRay;
			}

			points.push_back(lensRay.org);

			{
				double n1 = 1.0;
				for (int i = 0; i < m_lensSurface.size(); ++i) {
					auto& targetLens = m_lensSurface[i];

					std::optional<png::vec3> hitpoint;
					png::vec3 lensNormal;
					Intersect(png::vec3(0, 0, targetLens.z + targetLens.r), targetLens, lensRay, hitpoint, lensNormal);
					if (!hitpoint.has_value()) {
						// no ray
						generatedRay = std::nullopt;
						//points.push_back(lensRay.org + lensRay.dir);
						//Output(points, output);
						points.clear();
						return;
					}

					auto refractedDir = Refract(lensRay.dir, lensNormal, n1, targetLens.ior);

					if (!refractedDir.has_value()) {
						// no ray
						generatedRay = std::nullopt;
						//Output(points, output);
						points.clear();
						return;
					}
					lensRay.dir = refractedDir.value();
					lensRay.org = hitpoint.value();
					points.push_back(hitpoint.value());
					n1 = targetLens.ior;
				}
			}
			//auto t = (0.0 - lensRay.org.x) / lensRay.dir.x;
			auto t = 100.0;
			points.push_back(lensRay.org + t * lensRay.dir);
			Output(points, output);


			//Output(points, output);
			//return;

			m_fovx = m_fovy = 1.0;
			Ray generatedRayWorldSpace;
			generatedRayWorldSpace.org = m_camX * m_fovx * lensRay.org.x + m_camY * m_fovy * lensRay.org.y + m_data.cameraOrigin;
			generatedRayWorldSpace.dir = Normalize(m_camX * m_fovx * lensRay.dir.x + m_camY * m_fovy * lensRay.dir.y + m_camZ * lensRay.dir.z);
			std::vector<vec3> finalPoints;
			finalPoints.push_back(generatedRayWorldSpace.org);
			finalPoints.push_back(generatedRayWorldSpace.org + 100.0 * generatedRayWorldSpace.dir);
			//Output(finalPoints, output);
			generatedRay = generatedRayWorldSpace;
		}
	}

	void LensSystem::RayDir(double x, double y, double z, double spectrum, std::optional<Ray>& generatedRay, Random& rand, std::ofstream& output) {
		std::vector<png::vec3> points;
		Ray lensRay;
		{
			lensRay.org = vec3(x
				, y
				, m_lensSurface[0].z + z);
			std::optional<Ray> outRay;
			double xi_theta = std::acos(std::sqrt(rand.RandomGenerate()));
			double xi_phi = 2.0 * std::numbers::pi * rand.RandomGenerate();
			lensRay.dir = vec3(0, 0, 1);

		}

		points.push_back(lensRay.org);

		{
			double n1 = 1.0;
			for (int i = 0; i < m_lensSurface.size(); ++i) {
				auto& targetLens = m_lensSurface[i];

				std::optional<png::vec3> hitpoint;
				png::vec3 lensNormal;
				Intersect(png::vec3(0, 0, targetLens.z + targetLens.r), targetLens, lensRay, hitpoint, lensNormal);
				if (!hitpoint.has_value()) {
					// no ray
					generatedRay = std::nullopt;
					//points.push_back(lensRay.org + lensRay.dir);
					//Output(points, output);
					return;
				}

				auto refractedDir = Refract(lensRay.dir, lensNormal, n1, targetLens.ior);

				if (!refractedDir.has_value()) {
					// no ray
					generatedRay = std::nullopt;
					//Output(points, output);
					return;
				}
				lensRay.dir = refractedDir.value();
				lensRay.org = hitpoint.value();
				points.push_back(hitpoint.value());
				n1 = targetLens.ior;
			}
		}

		auto t = (0.0 - lensRay.org.x) / lensRay.dir.x;
		points.push_back(lensRay.org + t * lensRay.dir);

		//Output(points, output);
		//return;

		m_fovx = m_fovy = 1.0;
		Ray generatedRayWorldSpace;
		generatedRayWorldSpace.org = m_camX * m_fovx * lensRay.org.x + m_camY * m_fovy * lensRay.org.y + m_data.cameraOrigin;
		generatedRayWorldSpace.dir = Normalize(m_camX * m_fovx * lensRay.dir.x + m_camY * m_fovy * lensRay.dir.y + m_camZ * lensRay.dir.z);
		std::vector<vec3> finalPoints;
		finalPoints.push_back(generatedRayWorldSpace.org);
		finalPoints.push_back(generatedRayWorldSpace.org + 100.0 * generatedRayWorldSpace.dir);
		Output(finalPoints, output);
		return;
		generatedRay = generatedRayWorldSpace;
	}
}