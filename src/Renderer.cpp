#include "Renderer.h"

#include <numbers>
#include <random>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb_image_write.h>

namespace png {
	Renderer::Renderer(nlohmann::json& json)
		: width(json["width"])
		, height(json["height"])
		, sample(json["sample"])
		, image(std::vector<unsigned char>(width* height * 3))
	{}

	vec3 PathTracing(vec3 o, vec3 d) {

	}

	void Renderer::Render(Camera& cam, std::string fileName) {
		//dir
		auto direction = vec3::Normalize(cam.GetTarget() - cam.GetOrigin());
		auto l_camX = vec3::Normalize(vec3::Cross(direction, cam.GetUpVec()));
		auto l_camY = vec3::Cross(l_camX, direction);
		auto l_camZ = direction;
		//fov
		double fovx = std::cos((90.0f - cam.GetFOV() / 2.0f) / 180 * std::numbers::pi);
		double fovy = fovx * height / width;
		//random
		std::random_device rnd;
		std::mt19937 mt(rnd());

		for (int s = 0; s < sample; ++s) {
			for (int y = 0; y < height; ++y) {
				for (int x = 0; x < width; ++x) {
					double yd = 2.0 * y / (height - 1) - 1.0;
					double xd = 2.0 * x / (width - 1) - 1.0;

					vec3 dir = vec3::Normalize(l_camX * fovx * (2.0f * ((float)x + mt()) / width - 1.0f) +
						l_camY * fovy * (2.0f * ((float)y + mt()) / height - 1.0f) +
						l_camZ);
					auto cal = PathTracing(cam.GetOrigin(), dir) / sample;
					image[x * 3 + y * width * 3] = (unsigned char)(255.0 * cal.x);
					image[x * 3 + y * width * 3 + 1] = (unsigned char)(255.0 * cal.y);
					image[x * 3 + y * width * 3 + 2] = (unsigned char)(255.0 * cal.z);
				}
			}
		}
		stbi_write_jpg(fileName.c_str(), width, height, 0, (void*)0, 0);
	}
}