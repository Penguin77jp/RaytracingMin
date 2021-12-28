#include "Texture.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

namespace png {
	Texture::Texture(){}

	TextureSolid::TextureSolid(vec3 color)
		: m_color(color) {}

	vec3 TextureSolid::GetColor(double x, double y) const {
		return m_color;
	}

	TextureImage::TextureImage(std::string fileName, float level = 1.0f)
	: m_level(level){
		data = stbi_load(fileName.c_str(), &width, &height, &bits, 0);
	}

	vec3 TextureImage::GetColor(double x, double y) const {
		const int pixelX = (int)(width * x);
		const int pixelY = (int)(height * (1.0 - y));

		return vec3(
			(double)data[pixelX*3 + pixelY * width*3 + 0] / 255.0,
			(double)data[pixelX*3 + pixelY * width*3 + 1] / 255.0,
			(double)data[pixelX*3 + pixelY * width*3 + 2] / 255.0
		) * m_level;
	}

	TextureImage::~TextureImage() {
		stbi_image_free(data);
	}
}