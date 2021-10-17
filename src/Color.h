#pragma once

#include <cmath>

#include "Constant.h"
#include "Ray.h"


namespace png {

	// 390nm ~830nm
	// 830 - 390 = 440

	/*
	for wavelengthNM = 0:step:1000
	shape = 10;
	RGB.x = gaussianISO(wavelengthNM, 640, shape);
	RGB.y = gaussianISO(wavelengthNM, 555, shape);
	RGB.z = gaussianISO(wavelengthNM, 455, shape);
	index = wavelengthNM / step + 1;
	spectrumColorWavelength(index) = wavelengthNM;
	XYZ = xyzFromRgb(RGB);
	spectrumColorXYZ(index) = XYZ.x + XYZ.y + XYZ.z;
	end

	figure
	title('hoge')
	hold on
	xlim([390 830])
	line(spectrumColorWavelength, spectrumColorXYZ, 'Color', 'black')

	for wavelengthNM = 0:step:1000
	wavelengthA = angstromFromNanometer(wavelengthNM);
	XYZ = xbybzbFromWavelength(wavelengthA);
	index = wavelengthNM / step + 1;
	wavelength(index) = wavelengthNM;
	r(index) = XYZ.x;
	g(index) = XYZ.y;
	b(index) = XYZ.z;
	end

	figure
	hold on
	xlim([390 830])
	line(wavelength, r, 'Color', 'red')
	line(wavelength, g, 'Color', 'green')
	line(wavelength, b, 'Color', 'blue')

	colorXYZ.x = 0;
	colorXYZ.y = 0;
	colorXYZ.z = 0;
	for wavelengthNM = 0:step:1000
	xbybzb = xbybzbFromWavelength(angstromFromNanometer(wavelengthNM));
	colorXYZ.x = colorXYZ.x + xbybzb.x * spectrumColorXYZ(1 + wavelengthNM / step);
	colorXYZ.y = colorXYZ.y + xbybzb.y * spectrumColorXYZ(1 + wavelengthNM / step);
	colorXYZ.z = colorXYZ.z + xbybzb.z * spectrumColorXYZ(1 + wavelengthNM / step);
	end
	"colorXyz"
	colorXYZ
	"xyz2rgb"
	rgbFromXyz(colorXYZ)
	"xyz2rgb by Matlab"
	xyz2rgb([colorXYZ.x colorXYZ.y colorXYZ.z])

	*/

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
		return  alpha * exp(-(t * t) / 2);
	}

	vec3 xbybzbFromWavelength(double lambda) {
		vec3 xyz;
		xyz.x = gaussianAISO(lambda, 1.056, 599.8, 37.9, 31.0) + gaussianAISO(lambda, 0.362, 442.0, 16.0, 26.7) + gaussianAISO(lambda, -0.065, 501.1, 20.4, 26.2);
		xyz.y = gaussianAISO(lambda, 0.821, 568.8, 46.9, 40.5) + gaussianAISO(lambda, 0.286, 530.9, 16.3, 31.1);
		xyz.z = gaussianAISO(lambda, 1.217, 437.0, 11.8, 36.0) + gaussianAISO(lambda, 0.681, 459.0, 26.0, 13.8);
		return xyz;
	}

    // value 0 ~ 1.8
    double spectrumFromRGB(vec3 rgb) {
        return -1;
    }

	/*
	function string = xyzString(XYZ)
			string = "(" + XYZ.x + " , " + XYZ.y + " , " + XYZ.z + ")";
	*/


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


		/*
		function XYZ = rbgbbbFromWavelength(lambda)
			XYZ = rgbFromXyz(xbybzbFromWavelength(lambda));
		end
		*/
}