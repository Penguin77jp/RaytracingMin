#pragma once

#include <optional>
#include <vector>

#include "Camera.h"
#include "SettingData.h"
#include "Transparent.h"
#include "Random.h"

namespace png {
	class Py_LensSystem : public Camera {
	public:
		Py_LensSystem(const std::string fileName, const SettingData& data);
		void GenerateRay(int x, int y, int superSampleX, int superSampleY, double spectrum, Ray& rayIncomingSensor, std::optional<Ray>& generatedRay, Random& rand);
		void Test(vec3 samplePoint, vec3 sampleDir, double spectrum, Ray& rayIncomingSensor, std::optional<Ray>& generatedRay, Random& rand, std::vector<vec3>& out_line, const double z);
		void RayDir(const vec3 incomingImageSensor, double spectrum, std::optional<Ray>& generatedRay, Random& rand, std::vector<png::vec3>& output);

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