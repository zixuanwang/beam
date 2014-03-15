#include "MCLT.h"

MCLT::MCLT(int M) : m_M(M), m_coeff(sqrt(2.0 / (double)m_M) / 2.0), m_ha(new double[2 * m_M]), m_u(new double[m_M]), m_v(new double[m_M]){
	/// reference: A Modulated Complex Lapped Transform and its Applications to Audio Processing, Malvar, 99
	// construct u[n]
	int m2 = m_M * 2;
	for (int n = 0; n < m2; ++n){
		m_ha[n] = -sin((n + 0.5) * PI / m2);
	}
}

MCLT::~MCLT(){

}

void MCLT::compute(double* input, std::complex<double>* output){
	int m2 = m_M * 2;
	int m_half = m_M / 2;
	for (int n = 0; n < m_half; ++n){
		m_u[n + m_half] = input[m2 - 1 - n] * m_ha[m2 - 1 - n] - input[n + m_M] * m_ha[n + m_M];
		m_u[m_half - 1 - n] = input[m_M - 1 - n] * m_ha[n] + input[n] * m_ha[m_M - 1 - n];
	}
	for (int n = 0; n < m_half; ++n){
		m_v[n + m_half] = input[m2 - 1 - n] * m_ha[m2 - 1 - n] + input[n + m_M] * m_ha[n + m_M];
		m_v[m_half - 1 - n] = -input[m_M - 1 - n] * m_ha[n] + input[n] * m_ha[m_M - 1 - n];
	}
	// compute DCT and DST
	fftw_plan p;
	p = fftw_plan_r2r_1d(m_M, m_u.get(), m_u.get(), FFTW_REDFT11, FFTW_ESTIMATE);
	fftw_execute(p);
	p = fftw_plan_r2r_1d(m_M, m_v.get(), m_v.get(), FFTW_RODFT11, FFTW_ESTIMATE);
	fftw_execute(p);
	fftw_destroy_plan(p);
	// copy to output
	for (int n = 0; n < m_M; ++n){
		output[n] = std::complex<double>(m_coeff * m_u[n], m_coeff * m_v[n]);
	}
}