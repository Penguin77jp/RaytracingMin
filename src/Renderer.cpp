#include "Renderer.h"
#include "Color.h"
#include "Spectrum.h"

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
	{
		waveSetting.push_back(std::pair(0.25, std::pair(10.0, 100.0)));
		
		waveSetting.push_back(std::pair(0.5, std::pair(5, 50)));
		waveSetting.push_back(std::pair(0.5, std::pair(5, 50)));
		waveSetting.push_back(std::pair(0.5, std::pair(5, 50)));
		waveSetting.push_back(std::pair(0.5, std::pair(5, 50)));
		waveSetting.push_back(std::pair(0.5, std::pair(5, 50)));
		waveSetting.push_back(std::pair(0.5, std::pair(5, 50)));
		waveSetting.push_back(std::pair(0.5, std::pair(5, 50)));
		waveSetting.push_back(std::pair(0.5, std::pair(5, 50)));
		waveSetting.push_back(std::pair(0.5, std::pair(5, 50)));
		waveSetting.push_back(std::pair(0.5, std::pair(5, 50)));

		waveSetting.push_back(std::pair(0.75, std::pair(5, 50)));
		waveSetting.push_back(std::pair(0.75, std::pair(5, 50)));
		waveSetting.push_back(std::pair(0.75, std::pair(5, 50)));
		waveSetting.push_back(std::pair(0.75, std::pair(5, 50)));
	}

	void Renderer::AnimationUpdate(int frame, int allFrame, const clock_t benchmarkTime, Random& random) {
		float lerp = (float)frame / allFrame;
		data.AnimationUpdate(lerp);

		// camera
		{
			float _p = 2.0 * std::numbers::pi * lerp;
			float _r = 0.80;
			data.cameraOrigin = vec3(_r * std::cos(_p), data.cameraOrigin.y, _r * std::sin(_p));

			//aperture
			double apertureSpeed;
			double target = 0.0;
			if (frame == 0)
				data.cameraAperture = 0.03;
			//data.cameraAperture = 0.1;
			if (lerp < 1.0 / 2) {
				apertureSpeed = 0.0015;
				target = 0.0;
			}
			else {
				apertureSpeed = 0.001;
				target = 0.01;
			}
			auto& _ap = data.cameraAperture;
			/*
			if (std::abs(_ap - target) > 1e-3)
				_ap += (target - _ap) / std::abs(_ap - target) * apertureSpeed;
			*/
			if (target > _ap) {
				_ap += apertureSpeed;
			}
			else if (_ap > target) {
				_ap -= apertureSpeed;
			}
		}

		// wave
		{
			for (int i = 0; i < waveSetting.size();) {
				auto& _wave = waveSetting[i];
				if (_wave.first <= lerp) {
					data.wave->m_wave.gaussSet(random.RandomGenerate(), random.RandomGenerate(), _wave.second.first, _wave.second.second);
					//std::cout << "wave " << i << std::endl;
					waveSetting.erase(waveSetting.begin() + i);
				}
				else {
					i++;
				}
			}
		}
	}

	double hit_sphere(const vec3 center, double radius, const Ray r) {
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

	vec3 Pathtracing(const Ray ray, SettingData& data, Random& rand, const int depth) {
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
				auto hitPoint = ray.org + dis * ray.dir;
				if (rand.RandomGenerate() <= kd(obj->color(hitPoint))) {
					HitRecord nextRec;
					auto nextRay = obj->ScatteredRay(ray, nextRec, -1, rand);
					auto nextPathTracing = Pathtracing(nextRay, data, rand, depth + 1);
					return weight * (obj->color(vec3()) / kd(obj->color(hitPoint)) * nextPathTracing + obj->emission(hitPoint));
				}
				else {
					return weight * obj->emission(hitPoint);
				}
			}
		}
		if (data.sceneLight != nullptr) {
			return data.sceneLight->GetColor(ray.dir);
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
				if (tmp_dis < dis && tmp_dis > 0) {
					dis = tmp_dis;
					hitObject = i;
				}
			}
		}

		if (hitObject != -1) {
			auto& obj = data.object[hitObject];
			if (dis > 0) {
				auto hitpoint = ray.org + dis * ray.dir;
				//double  AlbedoSpectrumValue = color::spectrumValueFromRGB(obj->color(hitpoint), spectrum);
				double  AlbedoSpectrumValue = SampledSpectrumFromRGB(obj->color(hitpoint), SpectrumType::Reflectance, spectrum);
				if (rand.RandomGenerate() <= AlbedoSpectrumValue) {
					HitRecord nextRec;
					auto nextRay = obj->ScatteredRay(ray, nextRec, spectrum, rand);
					auto nextPathTracing = SpectrumPathtracing(nextRay, spectrum, data, rand, depth + 1);
					//return weight * (nextPathTracing + color::spectrumValueFromRGB(obj->emission(hitpoint), spectrum));
					return weight * (nextPathTracing + SampledSpectrumFromRGB(obj->emission(hitpoint), SpectrumType::Illuminant, spectrum));
				}
				else {
					//return weight * color::spectrumValueFromRGB(obj->emission(hitpoint), spectrum);
					return weight * SampledSpectrumFromRGB(obj->emission(hitpoint), SpectrumType::Illuminant, spectrum);
				}
			}
		}
		return 0;
	}

	vec3 RenderPathtracing(int x, int y, int sx, int sy, int samples, SettingData& data, Camera* cam, Random& rnd) {
		vec3 cal;
		Ray rayIncomingSensor;
		Ray ray;
		for (int s = 0; s < samples; ++s) {
			cam->GenerateRay(x, y, sx, sy, -1, rayIncomingSensor, ray, rnd);
			cal += Pathtracing(ray, data, rnd, 0) / samples;
		}
		return cal;
	}

	/*
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
	*/

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

	/*
	vec3 RenderSpectrumPathtracing_old(int x, int y, int sx, int sy, int samples, int spectrumSamples, Camera* cam, SettingData& data, Random& rnd) {
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
	*/

	vec3 RenderSpectrumPathtracing(int x, int y, int sx, int sy, int samples, int spectrumSamples, Camera* cam, SettingData& data, Random& rnd) {
		auto xyz_normalize = color::CIEXYZ::y_integral;

		vec3 xyz;
		for (int spec = 0; spec < spectrumSamples; ++spec) {
			double spectrum = rnd.RandomGenerate() * (color::MAX_WAVELENGTH - color::MIN_WAVELENGTH) + color::MIN_WAVELENGTH;
			Ray rayIncomingSensor;
			Ray ray;
			cam->GenerateRay(x, y, sx, sy, spectrum, rayIncomingSensor, ray, rnd);
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
			double spectrum = (color::MAX_WAVELENGTH - color::MIN_WAVELENGTH) * rnd.RandomGenerate() + color::MIN_WAVELENGTH;
			Ray rayIncomingSensor;
			Ray ray;
			cam->GenerateRay(x, y, sx, sy, spectrum, rayIncomingSensor, ray, rnd);
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
			auto max = 2.165515444135838851735798016306944191455841064453125; // means max value of x+y+z
			double xi_y = rand.RandomGenerate() * max;
			auto point = color::xbybzbFromWavelength(xi_x);
			double y = point.x + point.y + point.z;
			if (xi_y <= y) {
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
			cam->GenerateRay(x, y, sx, sy, spectrum, rayIncomingSensor, ray, rnd);
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
		const int renderingMode = data.renderType;
		double thickness = 1;
		double radius = 2;
		double focuce = 10;
		//double refractiveIndex = (radius + 0.5 * thickness) / 2 / focuce + 1;
		//Camera* camera = new NoLensCamera(data);
		ThinLensCamera thinCam = ThinLensCamera(data);
		thinCam.SetFocalPoint(vec3(0, 0, 0));
		thinCam.SetAperture(data.cameraAperture);
		Camera* camera = &thinCam;
		//Camera* camera = new PrototypeCamera(1, 2, TransparentMaterialType::HighVariance, data);
		//random
		png::Random rnd;

		std::clock_t start = std::clock();

#ifdef _DEBUG
#else
#pragma omp parallel for schedule(dynamic)
#endif
		for (int y = 0; y < data.height; ++y) {
			auto progress = (double)y / data.height;
			//std::cout << y << " / " << data.height << " [left " << leftMin << "min]" << std::endl;

			for (int x = 0; x < data.width; ++x) {
				vec3 accumulatedColor = vec3(0, 0, 0);
				for (int sx = 1; sx <= data.superSamples; ++sx) {
					for (int sy = 1; sy <= data.superSamples; ++sy) {
						Ray rayIncomingSendor;
						vec3 cal;
						if (renderingMode == 0) {
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
		//std::cout << ((double)(std::clock() - start) / CLOCKS_PER_SEC) << "sec" << std::endl;

		// save image
		auto resultImage = std::vector<unsigned char>(data.width * data.height * 3);
		for (int i = 0; i < resultImage.size(); ++i) {
			resultImage[i] = (unsigned char)(255.0 * std::min(image[i], 1.0));
		}

		stbi_write_jpg((fileName + ".jpg").c_str(), data.width, data.height, 3, resultImage.data(), 100);
		//stbi_write_bmp((fileName + ".bmp").c_str(), data.width, data.height, 3, resultImage.data());

	}
	}