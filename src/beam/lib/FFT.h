#ifndef FFT_H_
#define FFT_H_

#include <string.h>
#include "GlobalConfig.h"
#include "Utils.h"

namespace Beam{
	class FFT {
	public:
		/// constructor with the size.
		FFT();
		~FFT();
		/// compute FFT. make sure that output is allocated.
		/// the input must have FRAME_SIZE.
		/// the output must have FRAME_SIZE.
		void analyze(std::vector<float>& input, std::vector<std::complex<float> >& output);
		/// compute IFFT. make sure that output is allocated.
		/// the input must have FRAME_SIZE.
		/// the output must have FRAME_SIZE.
		/// the output has half frame delay.
		void synthesize(std::vector<std::complex<float> >& input, std::vector<float>& output);
		/// implementation of fft. copy from aecfft.c
		static void AecCcsFwdFFT(float * xin, float * xout, bool coeffOrder);
		/// implementation of ifft. copy from aecfft.c
		static void AecCcsInvFFT(float * xin, float * xout, bool coeffOrder);
	private:
		float m_ha[FRAME_SIZE * 2];
		float m_input[FRAME_SIZE * 2];
		float m_input_prev[FRAME_SIZE];
		float m_prev[FRAME_SIZE];
		float m_current[FRAME_SIZE * 2];
		// data structures to compute FFT.
		static float wr_15[15];
		static float wi_15[15];
		static void FwdFFT_base15(float * xin, float * xout);
		static void InvFFT_base15(float * xin, float * xout);
	};
}

#endif /* FFT_H_ */
