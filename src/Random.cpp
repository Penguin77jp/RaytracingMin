#include "Random.h"
#include <limits>

namespace png {

    Random::Random() {
        if (type == 0) {
            generator = new std::random_device();
        }else if (type == 1) {
            auto seed = std::random_device();
            generator = new std::mt19937(seed());
        }
    }

    Random::~Random() {
        if (type == 0) {
            delete (std::random_device*)(generator);
        }else if (type == 1) {
            delete (std::mt19937*)(generator);
        }
    }

    double Random::RandomGenerate() {
        unsigned  int randomInt;
        if (type == 0) {
            std::random_device* instance = (std::random_device*)(generator);
            randomInt = instance->operator()();
        }else if (type == 1) {
            std::mt19937* instance = (std::mt19937*)(generator);
            randomInt = instance->operator()();
        }
        return (double)(randomInt) / std::numeric_limits<unsigned  int>::max();
    }
}