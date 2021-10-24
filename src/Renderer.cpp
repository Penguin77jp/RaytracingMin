#include "Renderer.h"
#include "Color.h"

#include <numbers>
#include <random>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb_image_write.h>
#include <iostream>
#include <omp.h>
#include <ctime>

namespace png {
	using RandType = std::mt19937;
	Renderer::Renderer(SettingData& data)
		:image(std::vector<double>(data.width* data.height * 3))
		, data(data)
	{}

	float hit_sphere(const vec3 center, double radius, const Ray r) {
		vec3 oc = r.org - center;
		auto a = Dot(r.dir, r.dir);
		auto b = 2.0 * Dot(oc, r.dir);
		auto c = Dot(oc, oc) - radius * radius;
		auto discriminant = b * b - 4 * a * c;
		if (discriminant < 0) {
			return -1.0;
		}
		else {
			return (-b - sqrt(discriminant)) / (2.0 * a);
		}
	}

	double random(std::random_device& gene) {
		return (double)gene() / std::numeric_limits<unsigned int>::max();
	}

	vec3 clampColor(const vec3& ref, double min, float max) {
		return vec3(
			std::min<double>(std::max<double>(ref.x, min), max),
			std::min<double>(std::max<double>(ref.y, min), max),
			std::min<double>(std::max<double>(ref.z, min), max)
		);
	}

	constexpr bool ISDEBUG() {
#ifdef _DEBUG
		return true;
#else
		return false;
#endif
	}

	vec3 Pathtracing(const Ray ray, SettingData& data, Random& rand) {
		int hitObject = -1;
		double dis = std::numeric_limits<double>::max();
		{
			for (int i = 0; i < data.object.size(); ++i) {
				auto& obj = data.object[i];
				auto tmp_dis = obj->HitDistance(ray);
				//float tmp_dis = hit_sphere(obj.position, obj.size, ray);
				if (tmp_dis < dis && tmp_dis > 0) {
					dis = tmp_dis;
					hitObject = i;
				}
			}
		}

		if (hitObject != -1) {
			auto& obj = data.object[hitObject];
			if (dis > 0) {
				if (rand.RandomGenerate() <= obj->material->kd()) {
					HitRecord nextRec;
					auto nextRay = obj->ScatteredRay(ray, nextRec,  -1, rand);
					auto nextPathTracing = Pathtracing(nextRay, data, rand);
					return obj->material->color() / obj->material->kd() * nextPathTracing + obj->material->emission();
				}
				else {
					return obj->material->emission();
				}
			}
		}
		return vec3();
	}

	double SpectrumPathtracing(const Ray ray, const double spectrum, SettingData& data, Random& rand) {
		int hitObject = -1;
		double dis = std::numeric_limits<double>::max();
		{
			for (int i = 0; i < data.object.size(); ++i) {
				auto& obj = data.object[i];
				auto tmp_dis = obj->HitDistance(ray);
				//float tmp_dis = hit_sphere(obj.position, obj.size, ray);
				if (tmp_dis < dis && tmp_dis > 0) {
					dis = tmp_dis;
					hitObject = i;
				}
			}
		}

		if (hitObject != -1) {
			auto& obj = data.object[hitObject];
			if (dis > 0) {
				const double  spectrumValue = color::spectrumValueFromRGB(obj->material->color(), spectrum);
				if (rand.RandomGenerate() <= spectrumValue) {
					HitRecord nextRec;
					auto nextRay = obj->ScatteredRay(ray, nextRec, spectrum, rand);
					auto nextPathTracing = SpectrumPathtracing(nextRay, spectrum, data, rand);
					return nextPathTracing + color::spectrumValueFromRGB(obj->material->emission(), spectrum);
				}
				else {
					return color::spectrumValueFromRGB(obj->material->emission(), spectrum);
				}
			}
		}
		return 0;
	}

	vec3 RenderPathtracing(int samples, const Ray ray, SettingData& data, Random& rnd) {
		vec3 cal;
		for (int s = 0; s < samples; ++s) {
			cal += Pathtracing(ray, data, rnd) / samples;
		}
		return cal;
	}

	struct Spectrum {
		double spectrum, value;
		Spectrum(double spectrum, double value)
			: spectrum(spectrum)
			, value(value)
		{}
		Spectrum()
			: spectrum(0)
			, value(0)
		{}
	};
	bool operator< (const Spectrum& a, const Spectrum& b)  noexcept {
		return a.spectrum < b.spectrum;
	}

	vec3 RenderSpectrumPathtracing(int samples, int spectrumSamples, const Ray ray, SettingData& data, Random& rnd) {
		/*
		for (int wave = color::MIN_WAVELENGTH; wave <= color::MAX_WAVELENGTH; ++wave) {
			xyzNormalize += color::xbybzbFromWavelength(wave).y;
		}
		*/
		double xyzNormalize = 106.91612160775612;

		std::vector<Spectrum> spectrums(spectrumSamples);
		for (int spec = 0; spec < spectrumSamples; ++spec) {
			auto& cal = spectrums[spec];
			cal.spectrum = rnd.RandomGenerate() * (color::MAX_WAVELENGTH - color::MIN_WAVELENGTH) + color::MIN_WAVELENGTH;
			cal.value = 0;
			for (int s = 0; s < samples; ++s) {
				cal.value += SpectrumPathtracing(ray, cal.spectrum, data, rnd) / samples;
			}
		}

		std::sort(spectrums.begin(), spectrums.end());

		// spectrum to XYZ
		auto spectrumsTable = spectrums;
		spectrumsTable.insert(spectrumsTable.begin(), Spectrum(color::MIN_WAVELENGTH, 0));
		spectrumsTable.insert(spectrumsTable.end(), Spectrum(color::MAX_WAVELENGTH, 0));
		if (ISDEBUG()) {
			for (int i = 0; i < spectrumsTable.size(); ++i) {
//				std::cout << spectrumsTable[i].spectrum << "nm " << spectrumsTable[i].value << std::endl;
			}
		}


		vec3 XYZ;
		for (int index = 0; index < spectrumsTable.size() - 1; ++index) {
			double dx = spectrumsTable[index + 1].spectrum - spectrumsTable[index].spectrum;
			XYZ.x += (color::xbybzbFromWavelength(spectrumsTable[index + 1].spectrum).x * spectrumsTable[index + 1].value
				+ color::xbybzbFromWavelength(spectrumsTable[index].spectrum).x * spectrumsTable[index].value
				) * dx * 0.5 / xyzNormalize;
			XYZ.y += (color::xbybzbFromWavelength(spectrumsTable[index + 1].spectrum).y * spectrumsTable[index + 1].value
				+ color::xbybzbFromWavelength(spectrumsTable[index].spectrum).y * spectrumsTable[index].value
				) * dx * 0.5 / xyzNormalize;
			XYZ.z += (color::xbybzbFromWavelength(spectrumsTable[index + 1].spectrum).z * spectrumsTable[index + 1].value
				+ color::xbybzbFromWavelength(spectrumsTable[index].spectrum).z * spectrumsTable[index].value
				) * dx * 0.5 / xyzNormalize;
		}

		// XYZ to RGB
		auto RGB = color::rgbFromXyz(XYZ);

		return RGB;
	}

	void Renderer::Render(std::string fileName) {
		/*
		 *  Rendering Mode
		 *  0 : normal pathtracing
		 *  1: spectrum pathtracing
		 */

		const int renderingMode = 0;
		//dir
		auto direction = Normalize(data.camera.target - data.camera.origin);
		auto l_camX = -Normalize(Cross(direction, data.camera.upVec));
		auto l_camY = Cross(l_camX, direction);
		auto l_camZ = direction;
		//fov
		double fovx = data.camera.fov;
		double fovy = fovx * data.height / data.width;
        //random
        png::Random rnd;

		std::clock_t start = std::clock();
		for (int y = 0; y < data.height; ++y) {
			auto progress = (double)y / data.height;
			auto eclipseMin = (double)(std::clock() - start) / CLOCKS_PER_SEC / 60;
			auto leftMin = (int)(eclipseMin / progress - eclipseMin);
			std::cout << y << " / " << data.height << " [left " << leftMin << "min]" << std::endl;
#ifdef _DEBUG
#else
#pragma omp parallel for
#endif
			for (int x = 0; x < data.width; ++x) {
				vec3 accumulatedColor = vec3(0, 0, 0);
				for (int sx = 1; sx <= data.superSamples; ++sx) {
					for (int sy = 1; sy <= data.superSamples; ++sy) {
						const float rate = 1.0 / (1 + data.superSamples);
						vec3 dir = Normalize(
							l_camX * fovx * (2.0f * ((double)x + rate * sx) / data.width - 1.0f) +
							l_camY * fovy * (2.0f * ((double)y + rate * sy) / data.height - 1.0f) +
							l_camZ
						);

						vec3 cal;
						if (renderingMode == 0) {
							cal = RenderPathtracing(data.samples, Ray(data.camera.origin, dir), data, rnd);
						}
						else if (renderingMode == 1) {
							cal = RenderSpectrumPathtracing(data.samples, data.spectrumSamples, Ray(data.camera.origin, dir), data,
								rnd);
						}
						cal = clampColor(cal, 0, 1);
						accumulatedColor += cal / data.superSamples / data.superSamples;
					}
				}
				image[x * 3 + y * data.width * 3] += accumulatedColor.x;
				image[x * 3 + y * data.width * 3 + 1] += accumulatedColor.y;
				image[x * 3 + y * data.width * 3 + 2] += accumulatedColor.z;
			}
		}
		auto resultImage = std::vector<unsigned char>(data.width * data.height * 3);
		for (int i = 0; i < resultImage.size(); ++i) {
			resultImage[i] = (unsigned char)255 * std::min(image[i], 1.0);
		}

		stbi_write_jpg((fileName + ".jpg").c_str(), data.width, data.height, 3, resultImage.data(), 60);
		stbi_write_bmp((fileName + ".bmp").c_str(), data.width, data.height, 3, resultImage.data());

	}
}