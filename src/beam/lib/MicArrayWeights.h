#ifndef MICARRAYWEIGHTS_H_
#define MICARRAYWEIGHTS_H_

#include <complex>
#include <string>
#include <vector>
#include "Coords.h"

namespace Beam{
	class MicArrayWeights {
	public:
		std::string mic_array_name;
		int num_channels;
		int num_frequency_bins;
		int num_beams;
		float correlation;
		RCoords* beams;
		float* frequencies;
		std::complex<float>* weights;
		std::complex<float>* dd;
	};
}



#endif /* MICARRAYWEIGHTS_H_ */
