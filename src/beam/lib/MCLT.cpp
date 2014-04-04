#include "MCLT.h"

namespace Beam{
	MCLT::MCLT(int M) : m_M(M), m_coeff(sqrt(2.0 / (double)m_M) / 2.0), m_ha(new double[2 * m_M]), m_u(new double[m_M]), m_v(new double[m_M]){
		/// reference: A Modulated Complex Lapped Transform and its Applications to Audio Processing, Malvar, 99
		// construct u[n]
		int m2 = m_M * 2;
		for (int n = 0; n < m2; ++n){
			m_ha[n] = sin((n + 0.5) * PI / m2);
		}
	}

	MCLT::~MCLT(){

	}

	void MCLT::compute(const std::vector<float>& input, std::vector<std::complex<float> >& output){
		int m_half = m_M / 2;
		int m_three_halves = 3 * m_half;
		double t = 0.0;
		for (int n = 0; n < m_M; ++n){
			m_u[n] = m_ha[m_half + n] * input[m_three_halves - n - 1];
		}
		std::copy(m_u.get(), m_u.get() + m_M, m_v.get());
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
		fftw_plan p;
		p = fftw_plan_r2r_1d(m_M, m_u.get(), m_u.get(), FFTW_REDFT11, FFTW_ESTIMATE);
		fftw_execute(p);
		fftw_destroy_plan(p);
		p = fftw_plan_r2r_1d(m_M, m_v.get(), m_v.get(), FFTW_RODFT11, FFTW_ESTIMATE);
		fftw_execute(p);
		fftw_destroy_plan(p);
		// copy to output
		for (int n = 0; n < m_M; ++n){
			output[n] = std::complex<float>((float)(m_u[n] * m_coeff), (float)(m_v[n] * m_coeff));
		}
	}
}

