#include "FFT.h"

namespace Beam{
	FFT::FFT(int n) : m_n(n){
		for (int i = 0; i < m_n; ++i){
			m_ha[i] = sin((i + 0.5) * PI / n);
		}
	}
	FFT::~FFT(){
	
	}

	void FFT::compute(const std::vector<float>& input, std::vector<std::complex<float> >& output){
		fftw_complex in[FRAME_SIZE];
		fftw_complex out[FRAME_SIZE];
		for (int i = 0; i < FRAME_SIZE; ++i){
			in[i][0] = (double)input[i] * m_ha[i];
			in[i][1] = 0.0;
		}
		fftw_plan p;
		p = fftw_plan_dft_1d(FRAME_SIZE, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
		fftw_execute(p);
		fftw_destroy_plan(p);
		for (int i = 0; i < FRAME_SIZE; ++i){
			output[i] = std::complex<float>((float)out[i][0], (float)out[i][1]);
		}
	}
}

