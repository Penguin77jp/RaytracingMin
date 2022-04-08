#pragma once

#include "Ray.h"
#include "SceneObject.h"

#include <memory>

namespace png {
	struct HitRecord {
		vec3 point;
		double t;
		vec3 normal;
		//std::shared_ptr<Material> mat_ptr;
		bool front_face;
		double spectrum;
		int hitObjectIndex;

		inline void set_face_normal(const Ray& r, const vec3& outward_normal) {
			front_face = Dot(r.dir, outward_normal) < 0;
			normal = front_face ? outward_normal : -outward_normal;
		}
	};
}