#pragma once

#include <optional>
#include <vector>

#include "SettingData.h"
#include "Transparent.h"
#include "Random.h"

namespace png {
	struct SettingData;
	struct LensData {
		vec3 posi;
		double thickness;
		double radius_front, radius_rear;
		TransparentMaterialType traType;
	};

	class Camera {
	public:
		Camera(const SettingData& data);
		virtual void GenerateRay(int x, int y, int superSampleX, int superSampleY, double spectrum, Ray& rayIncomingSensor, std::optional<Ray>& generatedRay, Random& rand) = 0;

		//getter
		vec3 target() const;
		vec3 origin() const;
		double fov() const;
		int superSamples() const;

	public:
		const SettingData& m_data;
		vec3 m_target;
		vec3 m_origin;
		double m_fov;
		int m_superSamples;
	};

	class NoLensCamera : public Camera {
	public :
		NoLensCamera(SettingData& data);
		void GenerateRay(int x, int y, int superSampleX, int superSampleY, double spectrum, Ray& rayIncomingSensor, std::optional<Ray>& generatedRay, Random& rand);
	private:
		vec3 m_direction;
		vec3 m_camX;
		vec3 m_camY;
		vec3 m_camZ;
		double m_fovx;
		double m_fovy;
	};

	enum PrototypeLensType {
		Convex,
		Concave,
	};
	class LensCamera : public Camera {
	public:
		LensCamera(const std::vector<LensData>& lensData, const SettingData& data);
		LensCamera(const SettingData& data, const PrototypeLensType prototypeType);
		void GenerateRay(int x, int y, int superSampleX, int superSampleY, double spectrum, Ray& rayIncomingSensor, std::optional<Ray>& generatedRay, Random& rand);
	private:
		std::vector<LensData> m_lensData;
		vec3 m_aperturePosition;
		double m_apertureRadius;

		vec3 m_direction;
		vec3 m_camX;
		vec3 m_camY;
		vec3 m_camZ;
		double m_fovx;
		double m_fovy;
	};

	struct LensSurface {
		// pp : principal point
		double r, h, d, ior, z, f, pp;
	};

	class LensSystem : public Camera {
	public :
		LensSystem(const std::string fileName, const SettingData& data);
		void GenerateRay(int x, int y, int superSampleX, int superSampleY, double spectrum, Ray& rayIncomingSensor, std::optional<Ray>& generatedRay, Random& rand);
		void Test(int x, int y, double z, double spectrum, bool twoD, Ray& rayIncomingSensor, std::optional<Ray>& generatedRay, Random& rand, std::ofstream& output);
		void RayDir(double x, double y, double z, double spectrum, std::optional<Ray>& generatedRay, Random& rand, std::ofstream& output);

	public:
		std::vector<LensSurface> m_lensSurface;
		vec3 m_direction;
		vec3 m_camX;
		vec3 m_camY;
		vec3 m_camZ;
		double m_fovx;
		double m_fovy;
		double m_targetDistance;
	};
}