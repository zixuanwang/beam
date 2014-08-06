#ifndef GSCBEAMFORMER_H_
#define GSCBEAMFORMER_H_

#include "KinectConfig.h"
#include "SoundSourceLocalizer.h"

namespace Beam{
	class GSCBeamformer {
	public:
		GSCBeamformer();
		~GSCBeamformer();
		void compute(std::vector<std::complex<float> >* input, std::vector<std::complex<float> >& output, float angle, float confidence, double time, bool voice);
	private:
		std::vector<std::complex<float> > m_bm[MAX_MICROPHONES]; // store weights for blocking matrix
		std::vector<std::complex<float> > m_y[MAX_MICROPHONES]; // output of the blocking matrix
		std::vector<std::complex<float> > m_w[MAX_MICROPHONES]; // store weights for input canceller
	};
}

#endif /* GSCBEAMFORMER_H_ */
