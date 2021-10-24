#pragma once

#include <random>

namespace png {
    class Random{
    public :
        Random();
        ~Random();
        double RandomGenerate();

    private:
        void* generator;
        int type = 1;
    };

}