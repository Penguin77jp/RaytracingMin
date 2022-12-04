#pragma once

#include "ray.h"
#include <stb_image.h>
#include <cmath>
#include <numbers>
#include <vector>
#include <iostream>

namespace png {
    class SceneLight {
    public:
        float multiply;

    private:
        vec3* ambientLight;
        std::string imagePath;
        unsigned char* enviromentLight;
        int width, height, bits;

    public:
        SceneLight(float mul = 1.0f) 
            : multiply(mul)
            , ambientLight(nullptr)
            , enviromentLight(nullptr)
        {}
        SceneLight(vec3* col, float mul = 1.0f) 
            : multiply(mul)
            , ambientLight(col)
            , enviromentLight(nullptr)
        {}
        SceneLight(std::string imagePath, float mul = 1.0f) 
            : multiply(mul)
            , ambientLight(nullptr)
            , imagePath(imagePath)
            , enviromentLight(nullptr)
        {
            enviromentLight = stbi_load(imagePath.c_str(), &width, &height, &bits, 0);
        }

        vec3 GetColor(vec3 dir) {
            if (ambientLight != nullptr) {
                return *ambientLight * multiply;
            }
            else if (enviromentLight != nullptr) {
                const float length = std::sqrtf(dir.x * dir.x + dir.y * dir.y + dir.z * dir.z);
                double theta = std::acosf(dir.y / length);
                double phi = dir.x / std::abs(dir.x) * std::acosf(dir.z / std::sqrtf(dir.z * dir.z + dir.x * dir.x));
                double theta01 = theta / std::numbers::pi;
                double phi01 = phi / std::numbers::pi * 0.5f + 0.5f;
                if (std::isnan(phi01) || std::abs(phi01 - 1.0f) < 1.0e-10f) {
                    phi01 = 0;
                }
                if (std::isnan(theta01) || std::abs(theta01 - 1.0f) < 1.0e-10f) {
                    theta01 = 0;
                }
                int indX = phi01 * width;
                int indY = theta01 * height;
                int index = indY * width * 3 + indX * 3;
                index = std::min(index, 3*width*height);
                return vec3(enviromentLight[index], enviromentLight[index+1], enviromentLight[index+2]) / 255 * multiply;
            }
            else {
                return vec3{};
            }
        }

    };
}