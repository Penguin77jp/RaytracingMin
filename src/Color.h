#pragma once

#include <cmath>

#include "Constant.h"
#include "Ray.h"


namespace png {
    namespace color {
        // 390nm ~830nm
        // 830 - 390 = 440
        constexpr double MIN_WAVELENGTH = 390;
        constexpr double MAX_WAVELENGTH = 830;

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

        double gaussianISO(double x, double mu, double sigma);

        double gaussianAISO(double x, double alpha, double mu, double sigma1, double sigma2);

        vec3 xbybzbFromWavelength(double lambda);

        vec3 rgbFromXyz(vec3 XYZ);

        vec3 xyzFromRgb(vec3 RGB);

        double spectrumValueFromRGB(const vec3 rgb, const double wavelength);
    }
}