#include "FFT.h"

namespace Beam{
	FFT::FFT(){
		int two_frame_size = 2 * FRAME_SIZE;
		for (int i = 0; i < two_frame_size; ++i){
			m_ha[i] = sinf(i * (float)PI / (two_frame_size - 1));
		}
	}
	FFT::~FFT(){
	
	}

	void FFT::analyze(std::vector<float>& input, std::vector<std::complex<float> >& output){
		int two_frame_size = 2 * FRAME_SIZE;
		for (int i = 0; i < two_frame_size; ++i){
			m_input[i] = input[i] * m_ha[i];
		}
		fftwf_complex out[FRAME_SIZE + 1];
		fftwf_plan p = fftwf_plan_dft_r2c_1d(two_frame_size, m_input, out, FFTW_ESTIMATE);
		fftwf_execute(p);
		fftwf_destroy_plan(p);
		// only copy the first half.
		for (int i = 0; i <= FRAME_SIZE; ++i){
			output[i].real(out[i][0]);
			output[i].imag(out[i][1]);
		}
	}

	void FFT::synthesize(std::vector<std::complex<float> >& input, std::vector<float>& output){
		fftwf_complex in[2 * FRAME_SIZE];
		for (int i = 0; i <= FRAME_SIZE; ++i){
			in[i][0] = input[i].real();
			in[i][1] = input[i].imag();
		}
		fftwf_plan p = fftwf_plan_dft_c2r_1d(2 * FRAME_SIZE, in, &output[0], FFTW_ESTIMATE);
		fftwf_execute(p);
		fftwf_destroy_plan(p);
		
		memcpy(m_prev_output, &output[0], 2 * FRAME_SIZE * sizeof(float));
	}
}

