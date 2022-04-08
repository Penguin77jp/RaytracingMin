#pragma once

#include <random>

namespace png {
    class Random{
    public :
        Random();
        ~Random();
        double RandomGenerate();

    private:
        // 0 : std::random_device
        // 1 : std::mt19927
        int type = 1;
        void* generator;
    };

}