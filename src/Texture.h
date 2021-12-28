#pragma once

#include "Ray.h"

namespace png {

	class Texture {
	public:
		Texture();
		virtual vec3 GetColor(double x, double y) const = 0;
	};

	class TextureSolid : public Texture{
	public:
		TextureSolid(vec3 color);
		vec3 GetColor(double x, double y) const;
	private:
		vec3 m_color;
	};

	class TextureImage : public Texture {
	public:
		TextureImage(std::string fileName, float level);
		vec3 GetColor(double x, double y) const;
		~TextureImage();

	private:
		int width, height, bits;
		unsigned char* data;
		float m_level;
	};
}