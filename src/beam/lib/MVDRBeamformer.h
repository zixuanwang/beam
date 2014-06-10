#ifndef MVDRBEAMFORMER_H_
#define MVDRBEAMFORMER_H_

#include <armadillo>
#include "KinectConfig.h"
#include "SoundSourceLocalizer.h"

namespace Beam{
	class MVDRBeamformer {
	public:
		MVDRBeamformer();
		~MVDRBeamformer();
		void compute(std::vector<std::complex<float> >* input, std::vector<std::complex<float> >& output, float angle, float confidence, double time, bool voice = false);
	private:
		std::vector<arma::cx_fmat> m_nn;
	};
}

#endif /* MVDRBEAMFORMER_H_ */
