#include <array>
#include "Ray.h"

namespace png {
	static const int sampledLambdaStart = 400;
	static const int sampledLambdaEnd = 700;
	static const int nRGB2SpectSamples = 32;
	extern const double RGB2SpectLambda[32];
	extern const double RGBRefl2SpectWhite[32];
	extern const double RGBRefl2SpectCyan[32];
	extern const double RGBRefl2SpectMagenta[32];
	extern const double RGBRefl2SpectYellow[32];
	extern const double RGBRefl2SpectRed[32];
	extern const double RGBRefl2SpectGreen[32];
	extern const double RGBRefl2SpectBlue[32];
	extern const double RGBIllum2SpectWhite[32];
	extern const double RGBIllum2SpectCyan[32];
	extern const double RGBIllum2SpectMagenta[32];
	extern const double RGBIllum2SpectYellow[32];
	extern const double RGBIllum2SpectRed[32];
	extern const double RGBIllum2SpectGreen[32];
	extern const double RGBIllum2SpectBlue[32];

	enum class SpectrumType { Reflectance, Illuminant };

	struct Spectrum {
		std::array<double, 60> v;

		Spectrum() {
			for (int i = 0; i < 60; ++i) {
				v[i] = 0.0;
			}
		}

		Spectrum Clamp() const {
			Spectrum result;
			for (int i = 0; i < 60; ++i) {
				result.v[i] = std::max(this->v[i], 0.0);
			}
			return result;
		}
		Spectrum operator+=(const Spectrum& ref) {
			for (int i = 0; i < 60; ++i) {
				this->v[i] += ref.v[i];
			}
			return *this;
		}
		Spectrum operator*=(const double mul) {
			for (int i = 0; i < 60; ++i) {
				this->v[i] *= mul;
			}
			return *this;
		}

	};

	inline Spectrum operator*(const Spectrum& a, const double mul) {
		Spectrum result;
		for (int i = 0; i < 60; ++i) {
			result.v[i] = a.v[i] * mul;
		}
		return result;
	}
	inline Spectrum operator*(const double mul, const Spectrum& a) {
		return a * mul;
	}

	Spectrum SampledSpectrumFromRGB(const double rgb[3], SpectrumType type);
	double SampledSpectrumFromRGB(const png::vec3& rgb, SpectrumType type, double wavelength);

	inline double Lerp(double t, double v1, double v2) { return (1 - t) * v1 + t * v2; }

	inline double AverageSpectrumSamples(const double* lambda, const double* vals, int n,
		double lambdaStart, double lambdaEnd) {
		//for (int i = 0; i < n - 1; ++i) CHECK_GT(lambda[i + 1], lambda[i]);
		//CHECK_LT(lambdaStart, lambdaEnd);
		// Handle cases with out-of-bounds range or single sample only
		auto hoge = vals[2];
		if (lambdaEnd <= lambda[0]) return vals[0];
		if (lambdaStart >= lambda[n - 1]) return vals[n - 1];
		if (n == 1) return vals[0];
		double sum = 0;
		// Add contributions of constant segments before/after samples
		if (lambdaStart < lambda[0]) sum += vals[0] * (lambda[0] - lambdaStart);
		if (lambdaEnd > lambda[n - 1])
			sum += vals[n - 1] * (lambdaEnd - lambda[n - 1]);

		// Advance to first relevant wavelength segment
		int i = 0;
		while (lambdaStart > lambda[i + 1]) ++i;
		//CHECK_LT(i + 1, n);

		// Loop over wavelength sample segments and add contributions
		auto interp = [lambda, vals](double w, int i) {
			return Lerp((w - lambda[i]) / (lambda[i + 1] - lambda[i]), vals[i],
				vals[i + 1]);
		};
		for (; i + 1 < n && lambdaEnd >= lambda[i]; ++i) {
			double segLambdaStart = std::max(lambdaStart, lambda[i]);
			double segLambdaEnd = std::min(lambdaEnd, lambda[i + 1]);
			sum += 0.5 * (interp(segLambdaStart, i) + interp(segLambdaEnd, i)) *
				(segLambdaEnd - segLambdaStart);
		}
		return sum / (lambdaEnd - lambdaStart);
	}

	class SpectrumTable {
	private:
		SpectrumTable() {
			const int sampledLambdaStart = 400;
			const int sampledLambdaEnd = 700;
			const int nSpectralSamples = 60;
			const int nRGB2SpectSamples = 32;
			for (int i = 0; i < nSpectralSamples; ++i) {
				auto wl0 = Lerp(double(i) / double(nSpectralSamples),
					sampledLambdaStart, sampledLambdaEnd);
				auto wl1 = Lerp(double(i + 1) / double(nSpectralSamples),
					sampledLambdaStart, sampledLambdaEnd);
				rgbRefl2SpectWhite.v[i] =
					AverageSpectrumSamples(RGB2SpectLambda, RGBRefl2SpectWhite,
						nRGB2SpectSamples, wl0, wl1);
				rgbRefl2SpectCyan.v[i] =
					AverageSpectrumSamples(RGB2SpectLambda, RGBRefl2SpectCyan,
						nRGB2SpectSamples, wl0, wl1);
				rgbRefl2SpectMagenta.v[i] =
					AverageSpectrumSamples(RGB2SpectLambda, RGBRefl2SpectMagenta,
						nRGB2SpectSamples, wl0, wl1);
				rgbRefl2SpectYellow.v[i] =
					AverageSpectrumSamples(RGB2SpectLambda, RGBRefl2SpectYellow,
						nRGB2SpectSamples, wl0, wl1);
				rgbRefl2SpectRed.v[i] = AverageSpectrumSamples(
					RGB2SpectLambda, RGBRefl2SpectRed, nRGB2SpectSamples, wl0, wl1);
				rgbRefl2SpectGreen.v[i] =
					AverageSpectrumSamples(RGB2SpectLambda, RGBRefl2SpectGreen,
						nRGB2SpectSamples, wl0, wl1);
				rgbRefl2SpectBlue.v[i] =
					AverageSpectrumSamples(RGB2SpectLambda, RGBRefl2SpectBlue,
						nRGB2SpectSamples, wl0, wl1);

				rgbIllum2SpectWhite.v[i] =
					AverageSpectrumSamples(RGB2SpectLambda, RGBIllum2SpectWhite,
						nRGB2SpectSamples, wl0, wl1);
				rgbIllum2SpectCyan.v[i] =
					AverageSpectrumSamples(RGB2SpectLambda, RGBIllum2SpectCyan,
						nRGB2SpectSamples, wl0, wl1);
				rgbIllum2SpectMagenta.v[i] =
					AverageSpectrumSamples(RGB2SpectLambda, RGBIllum2SpectMagenta,
						nRGB2SpectSamples, wl0, wl1);
				rgbIllum2SpectYellow.v[i] =
					AverageSpectrumSamples(RGB2SpectLambda, RGBIllum2SpectYellow,
						nRGB2SpectSamples, wl0, wl1);
				rgbIllum2SpectRed.v[i] =
					AverageSpectrumSamples(RGB2SpectLambda, RGBIllum2SpectRed,
						nRGB2SpectSamples, wl0, wl1);
				rgbIllum2SpectGreen.v[i] =
					AverageSpectrumSamples(RGB2SpectLambda, RGBIllum2SpectGreen,
						nRGB2SpectSamples, wl0, wl1);
				rgbIllum2SpectBlue.v[i] =
					AverageSpectrumSamples(RGB2SpectLambda, RGBIllum2SpectBlue,
						nRGB2SpectSamples, wl0, wl1);
			}
		}
		~SpectrumTable() {}
	public:
		static SpectrumTable& getInstance();

		Spectrum rgbRefl2SpectWhite;
		Spectrum rgbRefl2SpectCyan;
		Spectrum rgbRefl2SpectBlue;
		Spectrum rgbRefl2SpectGreen;
		Spectrum rgbRefl2SpectMagenta;
		Spectrum rgbRefl2SpectRed;
		Spectrum rgbRefl2SpectYellow;
		Spectrum rgbIllum2SpectWhite;
		Spectrum rgbIllum2SpectCyan;
		Spectrum rgbIllum2SpectBlue;
		Spectrum rgbIllum2SpectGreen;
		Spectrum rgbIllum2SpectMagenta;
		Spectrum rgbIllum2SpectRed;
		Spectrum rgbIllum2SpectYellow;
	};
	inline SpectrumTable& SpectrumTable::getInstance() {
		static SpectrumTable obj;
		return obj;
	}

}