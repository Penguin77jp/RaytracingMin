#include "Camera.h"


#include "Ray.h"
#include "Random.h"
#include "color.h"
#include "Transparent.h"

//debug
#include <iostream>

namespace png {
	void kai(double a, double b, double c, std::optional<double>& t1, std::optional<double>& t2) {
		double b24ac = b * b - 4 * a * c;
		if (b24ac > 0) {
			t1 = (-b + sqrt(b24ac)) / (2 * a);
			t2 = (-b - sqrt(b24ac)) / (2 * a);
		}
		else if (b24ac == 0) {
			t1 = (-b) / (2 * a);
			t2 = std::nullopt;
		}
		else {
			t1 = std::nullopt;
			t2 = std::nullopt;
		}
	}
	std::optional<double> hitSphere(double thickness, double radius, Ray ray, double centerSphere) {
		double r = radius + thickness * 0.5;
		std::optional<double> t1, t2;
		kai(ray.dir.x * ray.dir.x + ray.dir.y * ray.dir.y + ray.dir.z * ray.dir.z,
			2 * (ray.org.x * ray.dir.x + ray.org.y * ray.dir.y + ray.org.z * ray.dir.z - centerSphere * ray.dir.z),
			ray.org.x * ray.org.x + ray.org.y * ray.org.y + ray.org.z * ray.org.z - 2 * centerSphere * ray.org.z + centerSphere * centerSphere - r * r
			, t1, t2);
		std::optional<double> t;
		if (t1 > 0 && t2 > 0) {
			if (t1 > t2) {
				t = t2;
			}
			else {
				t = t1;
			}
		}
		else if (t1 > 0 && t2 < 0) {
			t = t1;
		}
		else if (t1 < 0 && t2 > 0) {
			t = t2;
		}
		else {
			t = std::nullopt;
		}
		return t;
	}
	png::vec3 random_in_unit_sphere(png::Random& random) {
		while (true) {
			auto x = 2.0 * random.RandomGenerate() - 1.0;
			auto y = 2.0 * random.RandomGenerate() - 1.0;
			auto z = 2.0 * random.RandomGenerate() - 1.0;
			auto p = png::vec3(x, y, z);
			if (Magnitude(p) >= 1) continue;
			return p;
		}
	}
	bool lensInside(double thickness, double radius, double offset, png::vec3 point) {
		auto lensXMax = sqrt(std::pow(radius + thickness * 0.5, 2) - radius * radius);
		auto yMax = offset + thickness * 0.5;
		auto yMin = offset - thickness * 0.5;
		return (-lensXMax <= point.x && point.x <= lensXMax
			&& yMin <= point.y && point.y <= yMax);
	}
	std::optional<png::Ray> Lnes(double thickness, double radius, double refractiveIndex, png::Ray ray) {
		auto t = hitSphere(thickness, radius, ray, radius);
		if (!t.has_value())
			return std::nullopt;
		auto hitpoint = ray.org + t.value() * ray.dir;
		//% line
			//plot3([ray.org(1), hitpoint(1)], [ray.org(2), hitpoint(2)], [ray.org(3), hitpoint(3)])

		auto z = Normalize(hitpoint - vec3(0, 0, +radius));
		auto x = Normalize(Cross(z, vec3(1, 0, 0)));
		auto y = Normalize(Cross(z, x));

		/*
		% debug line
			% plot3([hitpoint(1), hitpoint(1) + z(1)], [hitpoint(2), hitpoint(2) + z(2)], [hitpoint(3), hitpoint(3) + z(3)], 'color', 'black')
			% plot3([hitpoint(1), hitpoint(1) + x(1)], [hitpoint(2), hitpoint(2) + x(2)], [hitpoint(3), hitpoint(3) + x(3)], 'color', 'black')
			% plot3([hitpoint(1), hitpoint(1) + y(1)], [hitpoint(2), hitpoint(2) + y(2)], [hitpoint(3), hitpoint(3) + y(3)], 'color', 'black')
			*/

		auto sinthetaIN = sqrt(pow2(Dot(x, ray.dir)) + pow2(Dot(y, ray.dir))); //% = cos(90 - x) = sinx
		auto dirLensSpaceIn = Normalize(vec3(Dot(x, ray.dir), Dot(y, ray.dir), Dot(z, ray.dir)));
		auto dirZ = refractiveIndex * (1.0 - pow2(sinthetaIN / refractiveIndex));
		auto dirLensSpace = Normalize(vec3(dirLensSpaceIn.x, dirLensSpaceIn.y, dirZ));

		ray.org = hitpoint;
		ray.dir = Normalize(x * dirLensSpace.x + y * dirLensSpace.y - z * dirLensSpace.z);

		//% inside
		t = hitSphere(thickness, radius, ray, -radius);
		if (!t.has_value())
			return std::nullopt;
		hitpoint = ray.org + t.value() * ray.dir;

		//% line
			//plot3([ray.org(1), hitpoint(1)], [ray.org(2), hitpoint(2)], [ray.org(3), hitpoint(3)])


		z = Normalize(hitpoint - vec3(0, 0, -radius));
		x = Normalize(Cross(z, vec3(1, 0, 0)));
		y = Normalize(Cross(z, x));

		/*
		% debug line
			% plot3([hitpoint(1), hitpoint(1) + z(1)], [hitpoint(2), hitpoint(2) + z(2)], [hitpoint(3), hitpoint(3) + z(3)], 'color', 'black')
			% plot3([hitpoint(1), hitpoint(1) + x(1)], [hitpoint(2), hitpoint(2) + x(2)], [hitpoint(3), hitpoint(3) + x(3)], 'color', 'black')
			% plot3([hitpoint(1), hitpoint(1) + y(1)], [hitpoint(2), hitpoint(2) + y(2)], [hitpoint(3), hitpoint(3) + y(3)], 'color', 'black')
			*/

		sinthetaIN = sqrt(pow2(Dot(x, ray.dir)) + pow2(Dot(y, ray.dir))); //% = cos(90 - x) = sinx
		dirLensSpaceIn = Normalize(vec3(Dot(x, ray.dir), Dot(y, ray.dir), Dot(z, ray.dir)));
		dirZ = (1 / refractiveIndex) * (1 - pow2((sinthetaIN * refractiveIndex)));
		dirLensSpace = Normalize(vec3(dirLensSpaceIn.x, dirLensSpaceIn.y, dirZ));

		//% outside
		ray.org = hitpoint;
		ray.dir = Normalize(x * dirLensSpace.x + y * dirLensSpace.y + z * dirLensSpace.z);

		//% line
			//t = 10;
		//plot3([ray.org(1), ray.org(1) + t * ray.dir(1)], [ray.org(2), ray.org(2) + t * ray.dir(2)], [ray.org(3), ray.org(3) + t * ray.dir(3)])
		return ray;
	}

	Camera::Camera(SettingData& data)
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
	void NoLensCamera::GenerateRay(int x, int y, int superSampleX, int superSampleY, double spectrum, Ray& rayIncomingSensor, Ray& generatedRay) {
		rayIncomingSensor = Ray();
		const float rate = 1.0 / (1 + superSamples());
		vec3 dir = Normalize(
			m_camX * m_fovx * (2.0f * ((double)x + rate * superSampleX) / m_data.width - 1.0f) +
			m_camY * m_fovy * (2.0f * ((double)y + rate * superSampleY) / m_data.height - 1.0f) +
			m_camZ
		);
		generatedRay = Ray(m_data.cameraOrigin, dir);
	}

	PrototypeCamera::PrototypeCamera(double thickness, double radius, TransparentMaterialType lensMaterialType, SettingData& data)
		: m_lensMaterialType(lensMaterialType)
		, m_thickness(thickness)
		, m_radius(radius)
		, Camera(data)
		, m_random(Random())
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

	void PrototypeCamera::GenerateRay(int x, int y, int superSampleX, int superSampleY, double spectrum, Ray& rayIncomingSensor, Ray& generatedRay) {
		const float rate = 1.0 / (1 + superSamples());
		Ray lensRay;
		lensRay.org = vec3((2.0 * ((double)x + rate * superSampleX) / m_data.width - 1.0)
			, 2.0 * ((double)y + rate * superSampleY) / m_data.height - 1.0, -2);
		std::optional<Ray> outRay;
		lensRay.dir = Normalize(vec3(0, 0, 1));
		//std::cout << refractiveIndex(m_lensMaterialType, spectrum) << std::endl;
		outRay = Lnes(m_thickness, m_radius, refractiveIndex(m_lensMaterialType, spectrum), lensRay);
		rayIncomingSensor = lensRay;


		//aperture
		// z = 10
		double apertureRadius = 0.12;
		double r = 0;
		{
			double t = (10.0 - outRay.value().org.z) / outRay.value().dir.z;
			auto point = outRay.value().org + t * outRay.value().dir;
			r = sqrt(point.x * point.x + point.y * point.y);
		}
		if (r >= apertureRadius) {
			//continue;
		}

		Ray ray;
		ray.org = origin() + m_camX * outRay.value().org.x + m_camY * outRay.value().org.y;
		//ray.org = origin();
		ray.dir = Normalize(m_camX * outRay.value().org.x + m_camY * outRay.value().org.y / m_data.width * m_data.height + m_camZ);

		generatedRay = ray;
	}
}