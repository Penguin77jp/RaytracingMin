#pragma once

#include "Ray.h"

namespace png {
	/*
	class Camera {
	public:
		Camera(vec3, vec3, vec3, float);
		vec3 origin, target, upVec;
		float fov;
	};
	*/
	class BiconvexLens  {
	public:
		BiconvexLens(float lensHeight, float thickness, float z);
		float lensHeight, thickness, z;

		// x = -1.0 <= x <= 1.0
		// y = -1.0 <= x <= 1.0
		Ray GetRay(double x, double y);
	};
}