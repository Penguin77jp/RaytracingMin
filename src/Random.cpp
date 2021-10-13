#include "Random.h"
#include <limits>

namespace png {

	double RandomGenerate(std::random_device& gene) {
		return (double)gene() / std::numeric_limits<unsigned int>::max();
	}

}