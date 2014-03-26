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
		unsigned long num_channels;
		unsigned long num_frequency_bins;
		unsigned long num_beams;
		float correlation;
		Coords* beams;
		float* frequencies;
		std::complex<float>* weights;
		std::complex<float>* dd;
	};
}



#endif /* MICARRAYWEIGHTS_H_ */
