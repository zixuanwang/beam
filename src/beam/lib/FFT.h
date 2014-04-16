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
		/// the output has half frame delay.
		void synthesize(std::vector<std::complex<float> >& input, std::vector<float>& output);
	private:
		float m_ha[FRAME_SIZE * 2];
		float m_input[FRAME_SIZE * 2];
		float m_prev_prev[FRAME_SIZE * 2];
		float m_prev[FRAME_SIZE * 2];
		float m_current[FRAME_SIZE * 2];
		// data structures to compute FFT.
		static float wr_15[15];
		static float wi_15[15];
		void FwdFFT_base15(float * xin, float * xout);
		void AecCcsFwdFFT(float * xin, float * xout, unsigned int FFTSize);
	};
}

#endif /* FFT_H_ */
