#ifndef DELAYSUMBEAMFORMER_H_
#define DELAYSUMBEAMFORMER_H_

#include "KinectConfig.h"
#include "SoundSourceLocalizer.h"

namespace Beam{
	class DelaySumBeamformer {
	public:
		DelaySumBeamformer();
		~DelaySumBeamformer();
		void compute(std::vector<std::complex<float> >* input, std::vector<std::complex<float> >& output, float angle, float confidence, double time);
	private:

	};
}

#endif /* DELAYSUMBEAMFORMER_H_ */
