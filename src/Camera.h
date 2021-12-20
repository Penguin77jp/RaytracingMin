#pragma once

#include <optional>

#include "SettingData.h"
#include "Transparent.h"

namespace png {
	struct SettingData;
	class Camera {
	public:
		Camera(SettingData& data);
		virtual void GenerateRay(int x, int y, int superSampleX, int superSampleY, double spectrum, Ray& rayIncomingSensor, Ray& generatedRay) = 0;

		//getter
		vec3 target() const;
		vec3 origin() const;
		double fov() const;
		int superSamples() const;

	protected:
		SettingData& m_data;
		vec3 m_target;
		vec3 m_origin;
		double m_fov;
		int m_superSamples;
	};

	class NoLensCamera : public Camera {
	public :
		NoLensCamera(SettingData& data);
		void GenerateRay(int x, int y, int superSampleX, int superSampleY, double spectrum, Ray& rayIncomingSensor, Ray& generatedRay);
	private:
		vec3 m_direction;
		vec3 m_camX;
		vec3 m_camY;
		vec3 m_camZ;
		double m_fovx;
		double m_fovy;
	};

	class PrototypeCamera : public Camera {
	public:
		PrototypeCamera(double thickness, double radius, TransparentMaterialType lensMaterialType, SettingData& data);
		void GenerateRay(int x, int y, int superSampleX, int superSampleY, double spectrum, Ray& rayIncomingSensor, Ray& generatedRay);
	private:
		TransparentMaterialType m_lensMaterialType;
		double m_thickness, m_radius;
		Random m_random;
		vec3 m_direction;
		vec3 m_camX;
		vec3 m_camY;
		vec3 m_camZ;
		double m_fovx;
		double m_fovy;
	};
}