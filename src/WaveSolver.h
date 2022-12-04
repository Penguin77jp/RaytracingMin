#pragma once

#include <vector>
#include <numbers>
#include <cmath>
#include "Ray.h"
#include "Random.h"

namespace {
	double distance(const png::vec3 a, const png::vec3 b) {
		return std::sqrt(std::pow(a.x - b.x, 2) + std::pow(a.y - b.y, 2));
	}
}

class WaveSolver {
public:
	WaveSolver() {
		_u_cur.resize(kWaveGrid * kWaveGrid);
		_u_new.resize(kWaveGrid * kWaveGrid);
		_u_pre.resize(kWaveGrid * kWaveGrid);
		this->reset();
	}
	void gaussSet(double _x, double _y, double sigma, double height) {
		double deltaX = double(kWaveWidth) / double(kWaveGrid);
		auto gauss = [](float x, float sigma) {
			return 1.0f / std::sqrt(2.0 * std::numbers::pi) * sigma * std::exp(-x * x / (2.0f * sigma * sigma));
		};
		for (int gx = 1; gx < kWaveGrid - 1; ++gx) {
			for (int gy = 1; gy < kWaveGrid - 1; ++gy) {
				int index = this->valueIndex(gx, gy);
				double norm = distance(png::vec3(gx, gy, 0) * deltaX, png::vec3(_x, _y, 0) * kWaveGrid * deltaX);
				_u_cur[index] += gauss(norm, sigma) * height;
			}
		}
		_u_pre = _u_new = _u_cur;
	}
	void reset() {
		gaussSet(0.5, 0.5, 10, 100);
	}
	void step() {
		double deltaX = double(kWaveWidth) / double(kWaveGrid);
		double c2 = 2.0f;
		double k = 0.03f;
		//float mul = deltaTime * deltaTime * c * c / (deltaX * deltaX);
#ifdef _DEBUG
#else
#pragma omp parallel for schedule(dynamic)
#endif
		for (int x = 1; x < kWaveGrid - 1; ++x) {
			for (int y = 1; y < kWaveGrid - 1; ++y) {
				

				int index = this->valueIndex(x, y);

				double uL = value(x - 1, y);
				double uR = value(x + 1, y);
				double uT = value(x, y - 1);
				double uB = value(x, y + 1);

				double u_pre = _u_pre[index];
				double u = _u_cur[index];

				//float val = u + u - u_pre + mul * (-4.0f * u + uL + uR + uT + uB) + damp;
				double val = c2 * png::pow2(deltaTime) / png::pow2(deltaX) *
					(-4.0 * u + uL + uR + uT + uB)
					- k * deltaTime * (u - u_pre) + 2.0 * u - u_pre;
				_u_new[index] = val;
			}
		}

		std::swap(_u_pre, _u_cur);
		std::swap(_u_cur, _u_new);
	}

	void step(int stepN) {
		for (int s = 0; s < stepN; ++s) {
			if (_rand.RandomGenerate() < 0.05) {
				//gaussSet(_rand.RandomGenerate(), _rand.RandomGenerate(), 35.0*_rand.RandomGenerate(), 250.0*std::pow(_rand.RandomGenerate(), 0.1));
			}
			step();
		}
	}

	enum {
		kWaveGrid = 400,
		//kWaveGrid = 1500,
		kWaveWidth = 500
	};
	double deltaTime = 1.0 / 10;

	int valueIndex(int x, int y) const {
		return y * kWaveGrid + x;
	}
	double value(int x, int y) const {
		if (x <= 0 || x >= widthN() || y <= 0 || y >= heightN()) {
			return 0;
		}
		return _u_cur[valueIndex(x, y)];
	}
	png::vec3 normal(int x, int y) const {
		auto p0 = value(x, y);
		auto pL = value(x - 1, y);
		auto pR = value(x + 1, y);
		auto pB = value(x, y - 1);
		auto pU = value(x, y + 1);
		const auto deltaX = width() / widthN();
		const auto deltaY = height() / heightN();
		auto crossed = Cross(
			png::vec3(0, (pR - pL) / (2.0 * deltaY),1),
			png::vec3(1, (pU - pB) / (2.0 * deltaX), 0)
		);
		crossed = Normalize(crossed);
		return crossed;
	}
	int widthN() const {
		return kWaveGrid;
	}
	int heightN() const {
		return kWaveGrid;
	}

	double width() const {
		return kWaveWidth;
	}
	double height() const {
		return kWaveWidth;
	}
	int time2step(double time) const {
		return time / deltaTime;
	}
private:
	std::vector<double> _u_pre;
	std::vector<double> _u_cur;
	std::vector<double> _u_new;
	png::Random _rand;
};
