#include "Color.h"


namespace png {
    namespace color {
        double gaussianISO(double x, double mu, double sigma) {
            double result = (1 / std::sqrt(2 * png::PI) / sigma) * exp(-((x - mu) * (x - mu)) / (2 * sigma * sigma));
            return result;
        }

        double gaussianAISO(double x, double alpha, double mu, double sigma1, double sigma2) {
            double sigma = 0;
            if (x < mu) {
                sigma = sigma1;
            }
            else {
                sigma = sigma2;
            }
            double t = (x - mu) / sigma;
            return alpha * exp(-(t * t) / 2);
        }

        vec3 xbybzbFromWavelength(double lambda) {
            vec3 xyz;
            xyz.x = gaussianAISO(lambda, 1.056, 599.8, 37.9, 31.0) + gaussianAISO(lambda, 0.362, 442.0, 16.0, 26.7) +
                gaussianAISO(lambda, -0.065, 501.1, 20.4, 26.2);
            xyz.y = gaussianAISO(lambda, 0.821, 568.8, 46.9, 40.5) + gaussianAISO(lambda, 0.286, 530.9, 16.3, 31.1);
            xyz.z = gaussianAISO(lambda, 1.217, 437.0, 11.8, 36.0) + gaussianAISO(lambda, 0.681, 459.0, 26.0, 13.8);
            return xyz;
        }

        vec3 rgbFromXyz(vec3 XYZ) {
            /*
                R = 3.2410X - 1.5374Y - 0.4986Z
                G = -0.9692X + 1.8760Y + 0.0416Z
                B = 0.0556X - 0.2040Y + 1.0507Z
            */
            vec3 RGB;
            RGB.x = 3.2410 * XYZ.x - 1.5374 * XYZ.y - 0.4986 * XYZ.z;
            RGB.y = -0.9692 * XYZ.x + 1.8760 * XYZ.y + 0.0416 * XYZ.z;
            RGB.z = 0.0556 * XYZ.x - 0.2040 * XYZ.y + 1.0507 * XYZ.z;
            return RGB;
        }

        vec3 xyzFromRgb(vec3 RGB) {
            /*
                X = 0.4124R + 0.3576G + 0.1805B
                Y = 0.2126R + 0.7152G + 0.0722B
                Z = 0.0193R + 0.1192G + 0.9505B
            */
            vec3 XYZ;
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


    }
}