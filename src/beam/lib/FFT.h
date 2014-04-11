#ifndef FFT_H_
#define FFT_H_

#include "fftw3.h"
#include "GlobalConfig.h"
#include "Utils.h"

namespace Beam{
	class FFT {
	public:
		/// constructor with the size.
		FFT();
		~FFT();
		/// compute FFT. make sure that output is allocated.
		/// the input must have 2 * FRAME_SIZE.
		/// the output must have FRAME_SIZE.
		void analyze(std::vector<float>& input, std::vector<std::complex<float> >& output);
		/// compute IFFT. make sure that output is allocated.
		/// the input must have FRAME_SIZE.
		/// the output must have 2 * FRAME_SIZE.
		void synthesize(std::vector<std::complex<float> >& input, std::vector<float>& output);
	private:
		float m_ha[FRAME_SIZE * 2];
		float m_input[FRAME_SIZE * 2];
		float m_prev_output[FRAME_SIZE * 2];
	};
}

#endif /* FFT_H_ */