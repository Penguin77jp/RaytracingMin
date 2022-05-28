#include <pybind11/pybind11.h>
#include <vector>
#include <pybind11/stl.h>
#include <optional>

#include "Ray.h"
#include "Camera_py.h"

using py_vec3 = std::vector<double>;

png::vec3 py_vec2vec(const py_vec3 x) {
	return png::vec3(x[0], x[1], x[2]);
}
py_vec3 vec2py_vec(const png::vec3 x) {
	return py_vec3{ x.x, x.y, x.z };
}


using py_line = std::vector<py_vec3>;

std::vector<py_line> cameraLensSample(const std::string settingDataPath, const std::string lensSystemPath, const py_vec3 samplePoint, const py_vec3 sampleDir, const double z) {
	png::Random rand;
	auto data = png::LoadData(settingDataPath);
	auto cam = png::Py_LensSystem(lensSystemPath, data.data);
	auto result = std::vector<py_line>();
	png::Ray incomingImageSensor;
	std::optional<png::Ray> generatedRay;
	std::vector<png::vec3> line;
	cam.Test(py_vec2vec(samplePoint), py_vec2vec(sampleDir), 0.0, incomingImageSensor, generatedRay, rand, line, z);
	//if (generatedRay.has_value()) {
		auto converted = py_line{};
		for (auto& p : line) {
			converted.push_back(vec2py_vec(p));
		}
		result.push_back(converted);
	//}
	return result;
}

std::vector<py_line> cameraLensSampleRayDir(const std::string settingDataPath, const std::string lensSystemPath, const py_vec3 samplePoint, const double samples) {
	png::Random rand;
	auto data = png::LoadData(settingDataPath);
	auto cam = png::Py_LensSystem(lensSystemPath, data.data);
	auto result = std::vector<py_line>();
	for (int i = 0; i < samples; ++i) {
		png::Ray incomingImageSensor;
		std::optional<png::Ray> generatedRay;
		std::vector<png::vec3> line;
		cam.RayDir(py_vec2vec(samplePoint), 0.0, generatedRay, rand, line);
		if (line.size() > 0) {
			auto converted = py_line{};
			for (auto& p : line) {
				converted.push_back(vec2py_vec(p));
			}
			result.push_back(converted);
		}
	}
	return result;
}

PYBIND11_MODULE(cppmodule, m) {
	m.def("cameraLensSample", &cameraLensSample);
	m.def("cameraLensSampleRayDir", &cameraLensSampleRayDir);
}
