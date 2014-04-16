#ifndef BEAMFORMER_H_
#define BEAMFORMER_H_

#include "KinectConfig.h"

namespace Beam{
#define F2RAISED23_INV (1.0f/8388608.0f)
	class Beamformer {
	public:
		Beamformer();
		~Beamformer();
		void init();
		void compute(std::vector<std::complex<float> >* input, std::vector<std::complex<float> >& output, float angle, double time);
		void ansi_bf_msr_process_quad_loop_fast(std::complex<float>* wo0, std::complex<float>* wo1, std::complex<float>* wo2, std::complex<float>* wo3, std::complex<float>& m0, std::complex<float>& m1, std::complex<float>& m2, std::complex<float>& m3, std::complex<float>& w0, std::complex<float>& w1, std::complex<float>& w2, std::complex<float>& w3, float nu, float mu);
	private:
		int m_beam;
		int m_first_bin;
		int m_last_bin;
		std::vector<std::complex<float> > m_pcm_weights[MAX_BEAMS][MAX_MICROPHONES];
	};
}

#endif /* BEAMFORMER_H_ */
