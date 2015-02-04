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

	// Fast forward MCLT implementation vis FFT
	// The windowing function must be strict sine window in this implementation
	// Reference: Fast Algorithm for the Modulated Complex Lapped Transform, Henrique S. Malvar, Technical Report, MSR-TR-2005-2
	// 
	// The output coefficient is arrange as below in case of DFT_COEFF_ORDER_AEC
	//      MCLT[0].re     --> AECFFTOut[0]
	//      MCLT[1].re     --> AECFFTOut[1]
	//        ...                ...
	//      MCLT[N/2-1].re --> AECFFTOut[N/2-1]
	//              0      --> AECFFTOut[N/2]
	//      MCLT[N/2-1].im --> AECFFTOut[N/2+1]
	//        ...                ...
	//      MCLT[2].im     --> AECFFTOut[N-2]
	//      MCLT[1].im     --> AECFFTOut[N-1]
	//
	// Note that MCLT[0].im is lost
	//
	// Buffer behavior: Input and output buffers can be same. If they are different, the input buffer data will not be changed. 
	void MCLT::AecCcsFwdMclt(float* pInput, float* pOutput, bool coeffOrder){
		int k, n, uFFTSize;
		float  r0, r1, i0, i1, ca, sa, uL, g;
		float  *up, tm, tp, cstep, sstep;
		float  * u;
		float  * y;

		float pfTempOut[FRAME_SIZE * 4];
		uFFTSize = FRAME_SIZE * 2;

		/* First compute FFT of input */
		// DFT_COEFF_ORDER_NRM has output with length of  FFTSize+2;
		FFT::AecCcsFwdFFT(pInput, pfTempOut, false);

		/* Get the size of the MCLT */
		n = uFFTSize / 2;
		u = pfTempOut;   // input of mapping at size of uFFTSize*2
		y = pOutput;    // output of mapping at size of uFFTSize

		/* Now apply DFT-to-MCLT mapping */
		uL = 0.5f / n;
		g = (float)(PI * (0.5f + uL));
		cstep = cosf(g);
		sstep = sinf(g);
		g = sqrtf(uL);
		ca = (float)cos(PI / 4.0);
		sa = -ca;
		r0 = u[0] * ca;
		i0 = u[0] * sa;
		up = &u[2];
		for (k = 0; k < n - 1; k++){
			r1 = *up++;
			i1 = *up++;
			tm = ca * cstep + sa * sstep;
			sa = sa * cstep - ca * sstep;
			ca = tm;
			tp = ca * r1 - sa * i1;
			i1 = sa * r1 + ca * i1;
			r1 = tp;
			y[k * 2] = g * (r1 - i0);
			y[k * 2 + 1] = g * (i1 + r0);
			r0 = r1; i0 = i1;
		}
		tm = ca * cstep + sa * sstep;
		sa = sa * cstep - ca * sstep;
		ca = tm;
		y[n * 2 - 2] = g * (ca * u[n * 2] - i0);
		y[n * 2 - 1] = g * (sa * u[n * 2] + r0);
		// re-oder the output for DFT_COEFF_ORDER_AEC
		if (coeffOrder)
		{
//			memcpy_s(pfTempOut, sizeof(float)*uFFTSize, pOutput, sizeof(float)*uFFTSize);
			memcpy(pfTempOut, pOutput, sizeof(float)*uFFTSize);
			pOutput[0] = pfTempOut[0];
			for (k = 1; k<uFFTSize / 2; k++)
			{
				pOutput[k] = pfTempOut[k * 2];
				pOutput[uFFTSize - k] = pfTempOut[k * 2 + 1];
			}
			pOutput[uFFTSize / 2] = 0;
		}
	}

	// Fast inverse MCLT implementation vis FFT
	// The windowing function must be strict sine window in this implementation
	// Reference: Fast Algorithm for the Modulated Complex Lapped Transform, Henrique S. Malvar, Technical Report, MSR-TR-2005-2
	// 
	// Buffer behavior: Input and output buffers can be same. If they are different, the input buffer data will not be changed. 
	void MCLT::AecCcsInvMclt(float * pInput, float * pOutput, bool coeffOrder)
	{
		int      k, n, j, uFFTSize;
		float   r1, i1, ca, sa, uL, g;
		float   tm, cstep, sstep;
		float   *y, *t;
		float   fltScale = 0.0f;

		float pfTempMCLTIn[FRAME_SIZE * 4];
		float pfTempFFTIn[FRAME_SIZE * 4];
		uFFTSize = FRAME_SIZE * 2;

		/* Get the size of the MCLT */
		n = uFFTSize / 2;

		// re-order input to normal order
		if (coeffOrder)
		{
			pfTempMCLTIn[0] = pInput[0];
			pfTempMCLTIn[1] = 0;
			for (k = 1; k < n; k++)
			{
				pfTempMCLTIn[k * 2] = pInput[k];
				pfTempMCLTIn[k * 2 + 1] = pInput[uFFTSize - k];
			}
			y = pfTempMCLTIn;
		}
		else{
			y = pInput;
		}

		t = pfTempFFTIn;

		/* Apply IMCLT-to-IDFT mapping */
		uL = 0.5f / n;
		g = (float)(PI * (0.5f + uL));
		cstep = cosf(g);
		sstep = sinf(-g);

		g = sqrtf(uL);
		ca = (float)cos(PI / 4.0);
		sa = ca;
		for (k = 1; k < n; k++)
		{
			r1 = y[k * 2 + 1] + y[k * 2 - 2];
			i1 = y[k * 2 - 1] - y[k * 2];
			tm = ca * cstep + sa * sstep;
			sa = sa * cstep - ca * sstep;
			ca = tm;
			t[k * 2] = ca * r1 - sa * i1;
			t[k * 2 + 1] = sa * r1 + ca * i1;
		}
		t[0] = sqrtf(2.f) * (y[0] + y[1]);
		t[1] = 0.f;
		t[n * 2] = -sqrtf(2.f) * (y[n * 2 - 2] + y[n * 2 - 1]);
		t[n * 2 + 1] = 0.f;
		k = n - 1;
		for (j = n + 1; j < 2 * n; j++)
		{
			t[j * 2] = t[k * 2];
			t[j * 2 + 1] = -t[k * 2 + 1];
			k--;
		}

		FFT::AecCcsInvFFT(pfTempFFTIn, pOutput, false);

		fltScale = 1.f / sqrtf(32.f * n);
		for (j = 0; j < uFFTSize; j++)
		{
			pOutput[j] *= fltScale;
		}
	}

	//void MCLT::analyze(std::vector<float>& input, std::vector<std::complex<float> >& output){
	//	int m_half = FRAME_SIZE / 2;
	//	int m_three_halves = 3 * m_half;
	//	float t = 0.f;
	//	for (int n = 0; n < FRAME_SIZE; ++n){
	//		m_u[n] = m_ha[m_half + n] * input[m_three_halves - n - 1];
	//	}
	//	std::copy(m_u, m_u + FRAME_SIZE, m_v);
	//	for (int n = 0; n < m_half; ++n){
	//		t = m_ha[m_half - 1 - n] * input[m_three_halves + n];
	//		m_u[n] += t;
	//		m_v[n] -= t;
	//	}
	//	for (int n = 0; n < m_half; ++n){
	//		t = m_ha[n] * input[n];
	//		m_u[n + m_half] -= t;
	//		m_v[n + m_half] += t;
	//	}
	//	// compute DCT and DST
	//	fftwf_plan p;
	//	p = fftwf_plan_r2r_1d(FRAME_SIZE, m_u, m_u, FFTW_REDFT11, FFTW_ESTIMATE);
	//	fftwf_execute(p);
	//	fftwf_destroy_plan(p);
	//	p = fftwf_plan_r2r_1d(FRAME_SIZE, m_v, m_v, FFTW_RODFT11, FFTW_ESTIMATE);
	//	fftwf_execute(p);
	//	fftwf_destroy_plan(p);
	//	// copy to output
	//	for (int n = 0; n < FRAME_SIZE; ++n){
	//		output[n].real(m_u[n] * m_coeff);
	//		output[n].imag(m_v[n] * m_coeff);
	//	}
	//}

	//void MCLT::synthesize(std::vector<std::complex<float> >& input, std::vector<float>& output){
	//	int m_half = FRAME_SIZE / 2;
	//	int m_three_halves = 3 * m_half;
	//	int m_twice = FRAME_SIZE * 2;
	//	for (int i = 0; i < FRAME_SIZE; ++i){
	//		m_u[i] = input[i].real();
	//		m_v[i] = input[i].imag();
	//	}
	//	// compute DCT and DST
	//	fftwf_plan p;
	//	p = fftwf_plan_r2r_1d(FRAME_SIZE, m_u, m_u, FFTW_REDFT11, FFTW_ESTIMATE);
	//	fftwf_execute(p);
	//	fftwf_destroy_plan(p);
	//	p = fftwf_plan_r2r_1d(FRAME_SIZE, m_v, m_v, FFTW_RODFT11, FFTW_ESTIMATE);
	//	fftwf_execute(p);
	//	fftwf_destroy_plan(p);
	//	for (int i = 0; i < m_half; ++i){
	//		m_current[i] = m_ha[i] * (m_v[m_half + i] - m_u[m_half + i]);
	//	}
	//	for (int i = m_half; i < m_three_halves; ++i){
	//		m_current[i] = m_ha[i] * (m_u[m_three_halves - i - 1] + m_v[m_three_halves - i - 1]);
	//	}
	//	for (int i = m_three_halves; i < m_twice; ++i){
	//		m_current[i] = m_ha[i] * (m_u[i - m_three_halves] - m_v[i - m_three_halves]);
	//	}
	//	for (int i = 0; i < m_twice; ++i){
	//		m_current[i] *= 0.5f * m_coeff;
	//	}
	//	// first half
	//	for (int i = 0; i < FRAME_SIZE; ++i){
	//		output[i] = m_prev_prev[i + FRAME_SIZE] + m_prev[i];
	//	}
	//	// second half
	//	for (int i = FRAME_SIZE; i < m_twice; ++i){
	//		output[i] = m_current[i - FRAME_SIZE] + m_prev[i];
	//	}
	//	// copy to prev prev
	//	memcpy(m_prev_prev, m_prev, m_twice * sizeof(float));
	//	// copy to prev
	//	memcpy(m_prev, m_current, m_twice * sizeof(float));
	//}
}

