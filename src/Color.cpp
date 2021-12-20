#include "Color.h"


namespace png {
	namespace color {
		double gaussianISO(double x, double mu, double sigma) {
			return (1 / std::sqrt(2 * png::PI) / sigma) * exp(-((x - mu) * (x - mu)) / (2 * sigma * sigma));
		}

		double gaussianAISO(double x, double mu, double sigma1, double sigma2) {
			double sigma = 0;
			if (x < mu) {
				sigma = sigma1;
			}
			else {
				sigma = sigma2;
			}
			double t = (x - mu) / sigma;
			return exp(-0.5 * (t * t));
		}

		vec3 xbybzbFromWavelength(double lambda) {
			vec3 xyz;
			xyz.x = 1.056 * gaussianAISO(lambda, 599.8, 37.9, 31.0) + 0.362 * gaussianAISO(lambda, 442.0, 16.0, 26.7) - 0.065 * gaussianAISO(lambda, 501.1, 20.4, 26.2);
			xyz.y = 0.821 * gaussianAISO(lambda, 568.8, 46.9, 40.5) + 0.286 * gaussianAISO(lambda, 530.9, 16.3, 31.1);
			xyz.z = 1.217 * gaussianAISO(lambda, 437.0, 11.8, 36.0) + 0.681 * gaussianAISO(lambda, 459.0, 26.0, 13.8);
			return xyz;
		}

		vec3 rgbFromXyz(vec3 XYZ) {
			vec3 RGB;
			/* CIE RGB
			RGB.x = 2.3655 * XYZ.x - 0.8971 * XYZ.y - 0.4683 * XYZ.z;
			RGB.y = -0.5151 * XYZ.x + 1.4264 * XYZ.y + 0.0887 * XYZ.z;
			RGB.z = 0.0052 * XYZ.x - 0.0144 * XYZ.y + 1.0089 * XYZ.z;
			*/
			// sRGB
			RGB.x = 3.2410 * XYZ.x - 1.5374 * XYZ.y - 0.4986 * XYZ.z;
			RGB.y = -0.9692 * XYZ.x + 1.8760 * XYZ.y + 0.0416 * XYZ.z;
			RGB.z = 0.0556 * XYZ.x - 0.2040 * XYZ.y + 1.0507 * XYZ.z;
			return RGB;
		}

		vec3 xyzFromRgb(vec3 RGB) {
			vec3 XYZ;
			/* CIE RGB
			XYZ.x = 0.4898 * RGB.x + 0.3101 * RGB.y + 0.2001 * RGB.z;
			XYZ.y = 0.1769 * RGB.x + 0.8124 * RGB.y + 0.0107 * RGB.z;
			XYZ.z = 0.0000 * RGB.x + 0.0100 * RGB.y + 0.9903 * RGB.z;
			*/

			//sRGB
			XYZ.x = 0.4124 * RGB.x + 0.3576 * RGB.y + 0.1805 * RGB.z;
			XYZ.y = 0.2126 * RGB.x + 0.7152 * RGB.y + 0.0722 * RGB.z;
			XYZ.z = 0.0193 * RGB.x + 0.1192 * RGB.y + 0.9505 * RGB.z;
			return XYZ;
		}

		double spectrumValueFromRGB(const vec3 rgb, const double wavelength) {
			const double r = 700;
			const double g = 546.1;
			const double b = 435.8;
			if (wavelength <= b) {
				return 0.0 + (wavelength - MIN_WAVELENGTH) / (b - MIN_WAVELENGTH) * (rgb.z - 0.0);
			}
			else if (wavelength <= g) {
				return rgb.z + (wavelength - b) / (g - b) * (rgb.y - rgb.z);
			}
			else if (wavelength <= r) {
				return rgb.y + (wavelength - g) / (r - g) * (rgb.x - rgb.y);
			}
			else {
				return rgb.x + (wavelength - r) / (MAX_WAVELENGTH - r) * (0.0 - rgb.x);
			}
		}

		namespace CIEXYZ {
			double gauseAISOCDF(double x, double alpha, double mu, double sigma1, double sigma2) {
				double sigma;
				if (x < mu) {
					sigma = sigma1;
				}
				else {
					sigma = sigma2;
				}
				return sqrt(png::PI / 2) * alpha * sigma * (erf(mu / (sqrt(2) * sigma)) + erf((x - mu) / (sqrt(2) * sigma)));
			}
		}
	}
}