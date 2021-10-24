#pragma once

#include <random>

namespace png {
    class Random{
    public :
        Random();
        ~Random();
        double RandomGenerate();
//        double RandomGenerate(std::random_device& gene);

    private:
        void* generator;
        int type = 1;
    };

	//using RandType = std::mt19937;
}