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
					auto nextRay = obj->ScatteredRay(ray, nextRec, -1, rand);
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

	double SpectrumPathtracing(const Ray ray, const double spectrum, SettingData& data, Random& rand, const int depth) {
		double weight = 1.0;
		{
			const double maxDepth = 1000;
			if (depth >= maxDepth) {
				const double expLambda = 1.0;
				double expVal = expLambda * exp(-expLambda * (depth - maxDepth));
				if (rand.RandomGenerate() > expVal) {
					return 0;
				}
				else {
					weight = 1.0 / std::min(expVal, 1.0);
				}
			}
		}
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
				double  AlbedoSpectrumValue = color::spectrumValueFromRGB(obj->material->color(), spectrum);
				if (rand.RandomGenerate() <= AlbedoSpectrumValue) {
					HitRecord nextRec;
					auto nextRay = obj->ScatteredRay(ray, nextRec, spectrum, rand);
					auto nextPathTracing = SpectrumPathtracing(nextRay, spectrum, data, rand, depth + 1);
					return weight * (nextPathTracing + color::spectrumValueFromRGB(obj->material->emission(), spectrum));
				}
				else {
					return weight * color::spectrumValueFromRGB(obj->material->emission(), spectrum);
				}
			}
		}
		return 0;
	}

	vec3 RenderPathtracing(int x, int y, int sx, int sy, int samples, SettingData& data, Camera* cam, Random& rnd) {
		vec3 cal;
		Ray rayIncomingSensor;
		Ray ray;
		cam->GenerateRay(x, y, sx, sy, -1, rayIncomingSensor, ray);
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

	double RandomSepctrum(Random& rnd) {
		namespace CIE = color::CIEXYZ;
		double random = rnd.RandomGenerate();
		double error = 0.1;
		double min = color::MIN_WAVELENGTH;
		double max = color::MAX_WAVELENGTH;
		/*
		while (true) {
			double mid = (min + max) * 0.5;
			double cal = (CIE::XCDF(mid) + CIE::YCDF(mid) + CIE::ZCDF(mid)) / CIE::XYZCDF_NORMALIZE;
			//std::cout << min << " , " << max << " -> " << std::abs(cal - random) << std::endl;
			if (std::abs(cal - random) <= error) {
				return mid;
			}
			else {
				if (cal > random) {
					max = mid;
				}
				else {
					min = mid;
				}
			}
		}
		*/
		return -1;
	}

	vec3 RenderSpectrumPathtracing_old(int x, int y, int sx, int sy, int samples, int spectrumSamples, Camera* cam, SettingData& data, Random& rnd) {
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
			//cal.spectrum = RandomSepctrum(rnd);
			Ray rayIncomingSensor;
			Ray ray;
			cam->GenerateRay(x, y, sx, sy, cal.spectrum, rayIncomingSensor, ray);
			cal.value = 0;
			for (int s = 0; s < samples; ++s) {
				//cal.value += SpectrumPathtracing(ray, cal.spectrum, data, rnd) / samples;
			}
			cal.value = 1.0;
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

	vec3 RenderSpectrumPathtracing(int x, int y, int sx, int sy, int samples, int spectrumSamples, Camera* cam, SettingData& data, Random& rnd) {
		auto xyz_normalize = color::CIEXYZ::y_integral;

		vec3 xyz;
		for (int spec = 0; spec < spectrumSamples; ++spec) {
			double spectrum = rnd.RandomGenerate() * (color::MAX_WAVELENGTH - color::MIN_WAVELENGTH) + color::MIN_WAVELENGTH;
			Ray rayIncomingSensor;
			Ray ray;
			cam->GenerateRay(x, y, sx, sy, spectrum, rayIncomingSensor, ray);
			for (int s = 0; s < samples; ++s) {
				xyz += (color::MAX_WAVELENGTH - color::MIN_WAVELENGTH) * color::xbybzbFromWavelength(spectrum) / samples / spectrumSamples / xyz_normalize;
				//xyz += (color::MAX_WAVELENGTH-color::MIN_WAVELENGTH) * SpectrumPathtracing(ray, spectrum, data, rnd)  * color::xbybzbFromWavelength(spectrum) / samples / spectrumSamples / xyz_normalize;
			}
		}

		auto RGB = color::rgbFromXyz(xyz);
		return RGB;
	}

	vec3 RenderSpectrumPathtracing_Normal(int x, int y, int sx, int sy, int samples, int spectrumSamples, Camera* cam, SettingData& data, Random& rnd) {
		vec3 xyz;
		for (int spec = 0; spec < spectrumSamples; ++spec) {
			double spectrum = (color::MAX_WAVELENGTH - color::MIN_WAVELENGTH)*rnd.RandomGenerate() + color::MIN_WAVELENGTH;
			Ray rayIncomingSensor;
			Ray ray;
			cam->GenerateRay(x, y, sx, sy, spectrum, rayIncomingSensor, ray);
			auto pathTracingRadiance = 0.0;
			for (int s = 0; s < samples; ++s) {
				pathTracingRadiance += SpectrumPathtracing(ray, spectrum, data, rnd, 0) / samples;
			}
			auto xbybzb = color::xbybzbFromWavelength(spectrum);
			xyz.x += pathTracingRadiance * xbybzb.x * 400.0 / color::CIEXYZ::y_integral / spectrumSamples;
			xyz.y += pathTracingRadiance * xbybzb.y * 400.0 / color::CIEXYZ::y_integral / spectrumSamples;
			xyz.z += pathTracingRadiance * xbybzb.z * 400.0 / color::CIEXYZ::y_integral / spectrumSamples;
		}

		auto RGB = color::rgbFromXyz(xyz);
		return RGB;
	}

	double SpectrumPoint_Importance(Random& rand) {
		while (true) {
			double xi_x = (color::MAX_WAVELENGTH - color::MIN_WAVELENGTH) * rand.RandomGenerate() + color::MIN_WAVELENGTH;
			auto max = 2.1655154441358; // means max value of x+y+z
			double xi_y = rand.RandomGenerate() * max;
			auto point = color::xbybzbFromWavelength(xi_x);
			double y = point.x + point.y + point.z;
			if (xi_y < y) {
				return xi_x;
			}
		}
	}

	vec3 RenderSpectrumPathtracing_ImportanceSampling(int x, int y, int sx, int sy, int samples, int spectrumSamples, Camera* cam, SettingData& data, Random& rnd) {
		vec3 xyz;
		for (int spec = 0; spec < spectrumSamples; ++spec) {
			double spectrum = SpectrumPoint_Importance(rnd);
			Ray rayIncomingSensor;
			Ray ray;
			cam->GenerateRay(x, y, sx, sy, spectrum, rayIncomingSensor, ray);
			auto pathTracingRadiance = 0.0;
			for (int s = 0; s < samples; ++s) {
				pathTracingRadiance += SpectrumPathtracing(ray, spectrum, data, rnd, 0) / samples;
			}
			auto xbybzb = color::xbybzbFromWavelength(spectrum);
			auto pdf = (xbybzb.x + xbybzb.y + xbybzb.z) / (color::CIEXYZ::x_integral + color::CIEXYZ::y_integral + color::CIEXYZ::z_integral);
			xyz.x += pathTracingRadiance * xbybzb.x / pdf / color::CIEXYZ::y_integral / spectrumSamples;
			xyz.y += pathTracingRadiance * xbybzb.y / pdf / color::CIEXYZ::y_integral / spectrumSamples;
			xyz.z += pathTracingRadiance * xbybzb.z / pdf / color::CIEXYZ::y_integral / spectrumSamples;
		}

		auto RGB = color::rgbFromXyz(xyz);
		return RGB;
	}

	void Renderer::Render(std::string fileName) {
		/*
		 *  Rendering Mode
		 *  0 : normal pathtracing
		 *  1 : spectrum pathtracing
		 */
		const int renderingMode = 0;
		double thickness = 1;
		double radius = 2;
		double focuce = 10;
		//double refractiveIndex = (radius + 0.5 * thickness) / 2 / focuce + 1;
		//Camera* camera = new NoLensCamera(data);
		Camera* camera = new PrototypeCamera(1, 2, TransparentMaterialType::HighVariance, data);
		//random
		png::Random rnd;

		std::clock_t start = std::clock();

#ifdef _DEBUG
#else
#pragma omp parallel for schedule(dynamic)
#endif
		for (int y = 0; y < data.height; ++y) {
			auto progress = (double)y / data.height;
			auto eclipseMin = (double)(std::clock() - start) / CLOCKS_PER_SEC / 60;
			auto leftMin = (int)(eclipseMin / progress - eclipseMin);
			std::cout << y << " / " << data.height << " [left " << leftMin << "min]" << std::endl;

			for (int x = 0; x < data.width; ++x) {
				vec3 accumulatedColor = vec3(0, 0, 0);
				for (int sx = 1; sx <= data.superSamples; ++sx) {
					for (int sy = 1; sy <= data.superSamples; ++sy) {
						Ray rayIncomingSendor;
						vec3 cal;
						if (renderingMode == 0) {
							cal = vec3(0, 1, 0);
							cal = RenderPathtracing(x, y, sx, sy, data.samples, data, camera, rnd);
						}
						else if (renderingMode == 1) {
							cal = RenderSpectrumPathtracing_ImportanceSampling(x, y, sx, sy, data.samples, data.spectrumSamples, camera, data, rnd);
						}
						//auto dot = Dot(rayIncomingSendor.dir, Normalize(vec3(0, 0, 1)));
						cal = clampColor(cal, 0, 1);
						accumulatedColor += cal / data.superSamples / data.superSamples;
					}
				}
				image[x * 3 + y * data.width * 3] += accumulatedColor.x;
				image[x * 3 + y * data.width * 3 + 1] += accumulatedColor.y;
				image[x * 3 + y * data.width * 3 + 2] += accumulatedColor.z;
			}
		}
		std::cout << ((double)(std::clock() - start) / CLOCKS_PER_SEC) << "sec" << std::endl;

		// save image
		auto resultImage = std::vector<unsigned char>(data.width * data.height * 3);
		for (int i = 0; i < resultImage.size(); ++i) {
			resultImage[i] = (unsigned char)255 * std::min(image[i], 1.0);
		}

		stbi_write_jpg((fileName + ".jpg").c_str(), data.width, data.height, 3, resultImage.data(), 60);
		stbi_write_bmp((fileName + ".bmp").c_str(), data.width, data.height, 3, resultImage.data());

	}
}