#pragma once

#include <vector>
#include <numbers>
#include <cmath>
#include "Ray.h"

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
	void reset() {
		auto gauss = [](float x, float sigma) {
			return 1.0f / std::sqrt(2.0 * std::numbers::pi) * sigma * std::exp(-x * x / (2.0f * sigma * sigma));
		};

		int cx1 = kWaveGrid / 3;
		int cy1 = kWaveGrid / 3;
		int cx2 = kWaveGrid * 2 / 3;
		int cy2 = kWaveGrid * 2 / 3;
		for (int x = 1; x < kWaveGrid - 1; ++x) {
			for (int y = 1; y < kWaveGrid - 1; ++y) {
				int index = this->valueIndex(x, y);
				float value = 0.0f;

				{
					float norm = distance(png::vec3(x, y, 0), png::vec3(cx1, cy1, 0));
					value += gauss(norm, 3.0f) * 20.0f;
				}
				{
					float norm = distance(png::vec3(x, y, 0), png::vec3(cx2, cy2, 0));
					value += gauss(norm, 3.0f) * 20.0f;
				}

				_u_cur[index] = value;
			}
		}

		_u_pre = _u_new = _u_cur;
	}
	void step() {
		float deltaX = float(kWaveWidth) / float(kWaveGrid);
		float c = 2.0f;
		float mul = deltaTime * deltaTime * c * c / (deltaX * deltaX);
		for (int x = 1; x < kWaveGrid - 1; ++x) {
			for (int y = 1; y < kWaveGrid - 1; ++y) {
				int index = this->valueIndex(x, y);

				float uL = _u_cur[this->valueIndex(x - 1, y)];
				float uR = _u_cur[this->valueIndex(x + 1, y)];
				float uT = _u_cur[this->valueIndex(x, y - 1)];
				float uB = _u_cur[this->valueIndex(x, y + 1)];

				float u_pre = _u_pre[index];
				float u = _u_cur[index];

				float k = 0.5f;
				float damp = -k * deltaTime * (u - u_pre);
				_u_new[index] = u + u - u_pre + mul * (-4.0f * u + uL + uR + uT + uB) + damp;
			}
		}

		std::swap(_u_pre, _u_cur);
		std::swap(_u_cur, _u_new);
	}

	void step(int stepN) {
		for (int s = 0; s < stepN; ++s) {
			step();
		}
	}

	enum {
		kWaveGrid = 300,
		kWaveWidth = 100
	};
	float deltaTime = 1.0 / 10;

	int valueIndex(int x, int y) const {
		return y * kWaveGrid + x;
	}
	float value(int x, int y) const {
		return _u_cur[valueIndex(x, y)];
	}
	png::vec3 normal(int x, int y) const {
		auto p0 = value(x, y);
		auto pL = value(x - 1, y);
		auto pR = value(x + 1, y);
		auto pB = value(x, y - 1);
		auto pU = value(x, y + 1);
		auto crossed = Cross(
			png::vec3(0, 1, (pU - pB) / 2),
			png::vec3(1, 0, (pR - pL) / 2)
		);
		return png::vec3(crossed.x, -crossed.z, crossed.y);
	}
	int widthN() const {
		return kWaveGrid;
	}
	int heightN() const {
		return kWaveGrid;
	}

	float width() const {
		return kWaveWidth;
	}
	float height() const {
		return kWaveWidth;
	}
	int time2step(float time) const {
		return time / deltaTime;
	}
private:
	std::vector<float> _u_pre;
	std::vector<float> _u_cur;
	std::vector<float> _u_new;
};
