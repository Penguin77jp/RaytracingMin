#include "Camera.h"
#include <iostream>

namespace png {
	BiconvexLens::BiconvexLens(float lensHeight, float thickness, float z)
		:lensHeight(lensHeight)
		, thickness(thickness)
		, z(z)
	{}

	Ray BiconvexLens::GetRay(double x, double y) {
		double circleCenter = (lensHeight * lensHeight - thickness * thickness / 4) / thickness;
		double theta = acos(circleCenter / (circleCenter + thickness / 2));
		double thetaX = (circleCenter + thickness / 2) * std::cos(theta * x);
		double thetaY = (circleCenter + thickness / 2) * std::cos(theta * y);
		vec3 Camera2LensDir = vec3(thetaX, thetaY, z);
		std::cout << Dot(-Normalize(Camera2LensDir),Normalize(vec3(thetaX, thetaY, -circleCenter))) << std::endl;
		return Ray();
	}
}