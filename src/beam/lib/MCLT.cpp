#include "MCLT.h"

namespace Beam{
	MCLT::MCLT() : m_coeff(sqrtf(2.f / (float)FRAME_SIZE) / 2.f){
		/// reference: A Modulated Complex Lapped Transform and its Applications to Audio Processing, Malvar, 99
		int two_frame_size = 2 * FRAME_SIZE;
		for (int i = 0; i < two_frame_size; ++i){
			m_ha[i] = sinf((0.5f + i) * (float)PI / two_frame_size);
			m_prev_prev[i] = 0.f;
			m_prev[i] = 0.f;
		}
	}

	MCLT::~MCLT(){

	}

	void MCLT::analyze(std::vector<float>& input, std::vector<std::complex<float> >& output){
		int m_half = FRAME_SIZE / 2;
		int m_three_halves = 3 * m_half;
		float t = 0.f;
		for (int n = 0; n < FRAME_SIZE; ++n){
			m_u[n] = m_ha[m_half + n] * input[m_three_halves - n - 1];
		}
		std::copy(m_u, m_u + FRAME_SIZE, m_v);
		for (int n = 0; n < m_half; ++n){
			t = m_ha[m_half - 1 - n] * input[m_three_halves + n];
			m_u[n] += t;
			m_v[n] -= t;
		}
		for (int n = 0; n < m_half; ++n){
			t = m_ha[n] * input[n];
			m_u[n + m_half] -= t;
			m_v[n + m_half] += t;
		}
		// compute DCT and DST
		fftwf_plan p;
		p = fftwf_plan_r2r_1d(FRAME_SIZE, m_u, m_u, FFTW_REDFT11, FFTW_ESTIMATE);
		fftwf_execute(p);
		fftwf_destroy_plan(p);
		p = fftwf_plan_r2r_1d(FRAME_SIZE, m_v, m_v, FFTW_RODFT11, FFTW_ESTIMATE);
		fftwf_execute(p);
		fftwf_destroy_plan(p);
		// copy to output
		for (int n = 0; n < FRAME_SIZE; ++n){
			output[n].real(m_u[n] * m_coeff);
			output[n].imag(m_v[n] * m_coeff);
		}
	}

	void MCLT::synthesize(std::vector<std::complex<float> >& input, std::vector<float>& output){
		int m_half = FRAME_SIZE / 2;
		int m_three_halves = 3 * m_half;
		int m_twice = FRAME_SIZE * 2;
		for (int i = 0; i < FRAME_SIZE; ++i){
			m_u[i] = input[i].real();
			m_v[i] = input[i].imag();
		}
		// compute DCT and DST
		fftwf_plan p;
		p = fftwf_plan_r2r_1d(FRAME_SIZE, m_u, m_u, FFTW_REDFT11, FFTW_ESTIMATE);
		fftwf_execute(p);
		fftwf_destroy_plan(p);
		p = fftwf_plan_r2r_1d(FRAME_SIZE, m_v, m_v, FFTW_RODFT11, FFTW_ESTIMATE);
		fftwf_execute(p);
		fftwf_destroy_plan(p);
		for (int i = 0; i < m_half; ++i){
			m_current[i] = m_ha[i] * (m_v[m_half + i] - m_u[m_half + i]);
		}
		for (int i = m_half; i < m_three_halves; ++i){
			m_current[i] = m_ha[i] * (m_u[m_three_halves - i - 1] + m_v[m_three_halves - i - 1]);
		}
		for (int i = m_three_halves; i < m_twice; ++i){
			m_current[i] = m_ha[i] * (m_u[i - m_three_halves] - m_v[i - m_three_halves]);
		}
		for (int i = 0; i < m_twice; ++i){
			m_current[i] *= 0.5f * m_coeff;
		}
		// first half
		for (int i = 0; i < FRAME_SIZE; ++i){
			output[i] = m_prev_prev[i + FRAME_SIZE] + m_prev[i];
		}
		// second half
		for (int i = FRAME_SIZE; i < m_twice; ++i){
			output[i] = m_current[i - FRAME_SIZE] + m_prev[i];
		}
		// copy to prev prev
		memcpy(m_prev_prev, m_prev, m_twice * sizeof(float));
		// copy to prev
		memcpy(m_prev, m_current, m_twice * sizeof(float));
	}
}

