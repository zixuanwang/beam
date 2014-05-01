#include "FFT.h"

namespace Beam{
	float FFT::wr_15[15] = {  // cos(2*pi*(0:14)/15)
		1.0f, 0.9135454576f, 0.6691306063f, 0.3090169943f, -0.1045284632f,
		-0.5f, -0.8090169943f, -0.9781476007f, -0.9781476007f, -0.8090169943f,
		-0.5f, -0.1045284632f, 0.3090169943f, 0.6691306063f, 0.9135454576f };
	float FFT::wi_15[15] = {  // -sin(2*pi*(0:14)/15)
		0.0f, -0.4067366430f, -0.7431448254f, -0.9510565162f, -0.9945218953f,
		-0.8660254037f, -0.5877852522f, -0.2079116908f, 0.2079116908f, 0.5877852522f,
		0.8660254037f, 0.9945218953f, 0.9510565162f, 0.7431448254f, 0.4067366430f };


	FFT::FFT(){
		int two_frame_size = 2 * FRAME_SIZE;
		for (int i = 0; i < two_frame_size; ++i){
			m_ha[i] = sinf(i * (float)PI / two_frame_size);
			m_prev_prev[i] = 0.f;
			m_prev[i] = 0.f;
		}
	}
	FFT::~FFT(){
	
	}

	void FFT::analyze(std::vector<float>& input, std::vector<std::complex<float> >& output){
		const int two_frame_size = 2 * FRAME_SIZE;
		float out[two_frame_size];
		for (int i = 0; i < two_frame_size; ++i){
			m_input[i] = input[i] * m_ha[i];
		}
		AecCcsFwdFFT(m_input, out, (unsigned int)two_frame_size);
		// only copy the first half
		for (int i = 0; i < FRAME_SIZE; ++i){
			output[i].real(out[i]);
		}
		for (int i = FRAME_SIZE + 1; i < two_frame_size; ++i){
			output[two_frame_size - i].imag(out[i]);
		}
		output[0].imag(0.f);
	}

	void FFT::synthesize(std::vector<std::complex<float> >& input, std::vector<float>& output){
		const int two_frame_size = 2 * FRAME_SIZE;
		float in[two_frame_size] = { 0.f };
		for (int i = 0; i < FRAME_SIZE; ++i){
			in[i] = input[i].real();
		}
		for (int i = FRAME_SIZE + 1; i < two_frame_size; ++i){
			in[i] = input[two_frame_size - i].imag();
		}
		AecCcsInvFFT(in, m_current, (unsigned int)two_frame_size);
		for (int i = 0; i < two_frame_size; ++i){
			m_current[i] /= (float)two_frame_size;
		}
		for (int i = 0; i < two_frame_size; ++i){
			m_current[i] *= m_ha[i];
		}
		// first half
		for (int i = 0; i < FRAME_SIZE; ++i){
			output[i] = m_prev_prev[i + FRAME_SIZE] + m_prev[i];
		}
		// second half
		for (int i = FRAME_SIZE; i < two_frame_size; ++i){
			output[i] = m_current[i - FRAME_SIZE] + m_prev[i];
		}
		// copy to prev prev
		memcpy(m_prev_prev, m_prev, two_frame_size * sizeof(float));
		// copy to prev
		memcpy(m_prev, m_current, two_frame_size * sizeof(float));
	}

	void FFT::FwdFFT_base15(float * xin, float * xout){
		int i, j;
		xout[0] = xin[0] + xin[1] + xin[2] + xin[3] + xin[4] + xin[5] + xin[6] + xin[7] +
			xin[8] + xin[9] + xin[10] + xin[11] + xin[12] + xin[13] + xin[14];
		for (i = 1; i<8; i++){
			xout[i] = xin[0];
			xout[15 - i] = 0;
			for (j = 1; j<15; j++){
				xout[i] += xin[j] * wr_15[(i*j) % 15];
				xout[15 - i] += xin[j] * wi_15[(i*j) % 15];
			}
		}
	}

	// Inverse base transform for FFT size of 15
	void FFT::InvFFT_base15(float * xin, float * xout)
	{
		int i, j;

		xout[0] = xin[0] + (xin[1] + xin[2] + xin[3] + xin[4] + xin[5] + xin[6] + xin[7]) * 2;
		for (i = 1; i<15; i++)
		{
			xout[i] = xin[0];
			for (j = 1; j<8; j++)
			{
				xout[i] += xin[j] * wr_15[(i*j) % 15] * 2;
				xout[i] += xin[15 - j] * wi_15[(i*j) % 15] * 2;
			}
		}
	}

	/***************************************************************************
	Function name: AecCcsFwdFFT

	Decription:
	This function supports two output coefficient orders - DFT_COEFF_ORDER_AEC and DFT_COEFF_ORDER_NRM.
	For DFT_COEFF_ORDER_AEC, pfltOutData is FFTSize long, otherwise it is FFTSize+2 long.

	pFltSinTable points to a FFTSize/4+1 long sin table which is necessary for ANSI C FFT function AecFwdFFT.
	IntelIPP function does not use it. Application is responsible for providing corrent sintable to this function.
	Please refer to filtbank.c for example.

	If FFTSize if power of 2 and DFT_COEFF_ORDER_AEC is specified, this function use in-place computation;
	otherwise a temp buffer will be used internally, but you can still use same buffer for input and
	output sequences. Please note that if DFT_COEFF_ORDER_NRM is specified, input length and output length are
	different. Make sure your buffer is long enough if you want to use in-place computation. If input and output
	buffers are different, the input buffer data will not be changed.

	Input:  pXformInfo:  AEC transform structure information
	xin:  Pointer of real time sequence in normal order x[0], x[1], ... , x[FFTSize-1];
	sin_table:  sin table containing values of sin(2*pi*(0:FFTSize/4)/FFTSize). The size
	of table is FFTSize/4+1;
	FFTSize:  length of time sequence, FFTSize must 2^n, 5*2^n, or 15*2^n
	fCoeffOrder: DFT_COEFF_ORDER_AEC or DFT_COEFF_ORDER_NRM

	Output: xout:  Pointer of complex frequency sequence, only  0 - FFTSize/2 complex samples
	are saved as follow. This order is also used for internal calculation
	in AEC. The orignal input is replaced by output.
	X[0], Xr[1], Xr[2], ..., Xr[FFTSize/2-1], X[FFTSize/2], Xi[FFTSize/2-1], ..., Xi[2], Xi[1];

	Implementation: Based on decimation-in-time radix-2 algorithm. First the FFT
	for size 4 or 5 is calcuated directly. Then the butterfly
	calculation is performed progressively.

	Reference:  Sorensen, Jones, at el, "Real Valued Fast Fourier Transform
	Algorithms," IEEE Trans. Acoust., Speech, Signal Processing,
	Vol ASSP-36, 1987.

	Qin Li, Feb 18, 2005
	***************************************************************************/
	void FFT::AecCcsFwdFFT(float * xin, float * xout, unsigned int FFTSize){
		// sin and cos table for base 5 FFT
		const float wr1 = 0.309016994374947f;  // cos(2pi/5)
		const float wr2 = -0.809016994374947f;  // cos(4pi/5)
		const float wr3 = -0.809016994374947f;  // cos(6pi/5)
		const float wr4 = 0.309016994374947f;  // cos(8pi/5)
		const float wi1 = -0.951056516295154f;  // -sin(2pi/5) 
		const float wi2 = -0.587785252292473f;  // -sin(4pi/5)
		const float wi3 = 0.587785252292473f;  // -sin(6pi/5)
		const float wi4 = 0.951056516295154f;  // -sin(4pi/5)

		float * tempbuf = new float[2 * FFTSize];
		float * sin_tab = new float[FFTSize / 4 + 1];
		for (unsigned int i = 0; i <= FFTSize / 4; i++) {
			sin_tab[i] = (float)sinf(2.0f * (float)PI * i / FFTSize);
		}
		unsigned int base, step, N_base;
		unsigned int i, j, k;   // loop indices
		float * cos_tab, *x;

		if (xin != xout) {
			memcpy(xout, xin, FFTSize*sizeof(float));
		}
		x = xout;

		cos_tab = sin_tab + FFTSize / 4;

		if ((FFTSize & (-(int)FFTSize)) == FFTSize) // detect if FFTSize is power of 2
		{
			base = 4;
		}
		else{
			if (FFTSize % 15)
				base = 5;
			else
				base = 15;
		}
		N_base = FFTSize / base;  //number of base transforms

		// bit reversal
		if (base == 4)  // FFTSize is power of 2
		{
			for (i = 0, j = 0; i < FFTSize; i++)
			{
				if (j>i)
				{
					float temp;
					temp = x[j];
					x[j] = x[i];
					x[i] = temp;
				}
				k = FFTSize / 2;
				while (k >= 2 && j >= k)
				{
					j -= k;
					k >>= 1;
				}
				j += k;
			}
		}
		else{
			if (base == 5)
			{
				/***********************************************************************
				cannot do in-place indexing.The basic idea is to do bit inverse
				for lower bits corresponding to N_base. The higest bits corresponding
				to the prime factor of 5 are not reversed. Since this is not pair-wise
				swap, it cannot be done in-place

				For example, for fft of 20 points
				normal order              bit-reversed order
				0   (00000)                 0   (00000)
				1   (00001)                 4   (00100)
				2   (00010)                 8   (01000)
				3   (00011)                 12  (01100)
				4   (00100)                 16  (10000)
				5   (00101)                 2   (00010)
				6   (00110)                 6   (00110)
				7   (00111)                 10  (01010)
				8   (01000)                 14  (01110)
				9   (01001)                 18  (10010)
				10  (01010)                 1   (00001)
				11  (01011)                 5   (00101)
				12  (01100)                 9   (01001)
				13  (01101)                 13  (01101)
				14  (01110)                 17  (10001)
				15  (01111)                 3   (00011)
				16  (10000)                 7   (00111)
				17  (10001)                 11  (01011)
				18  (10010)                 15  (01111)
				19  (10011)                 19  (10011)
				***********************************************************************/

				for (i = 0, j = 0; i < N_base && j < N_base; i++)
				{
					if (j == i)  //just copy, no bit reverse
					{
						unsigned int utemp;
						utemp = i*base;
						tempbuf[utemp] = x[i];
						tempbuf[utemp + 1] = x[N_base + i];
						tempbuf[utemp + 2] = x[N_base * 2 + i];
						tempbuf[utemp + 3] = x[N_base * 3 + i];
						tempbuf[utemp + 4] = x[N_base * 4 + i];
					}
					if (j > i)  // bit reverse and copy
					{
						unsigned int utemp;
						utemp = j*base;
						tempbuf[utemp] = x[i];
						tempbuf[utemp + 1] = x[N_base + i];
						tempbuf[utemp + 2] = x[N_base * 2 + i];
						tempbuf[utemp + 3] = x[N_base * 3 + i];
						tempbuf[utemp + 4] = x[N_base * 4 + i];

						utemp = i*base;
						tempbuf[utemp] = x[j];
						tempbuf[utemp + 1] = x[N_base + j];
						tempbuf[utemp + 2] = x[N_base * 2 + j];
						tempbuf[utemp + 3] = x[N_base * 3 + j];
						tempbuf[utemp + 4] = x[N_base * 4 + j];
					}
					k = N_base / 2;
					while (k >= 2 && j >= k)
					{
						j -= k;
						k >>= 1;
					}
					j += k;
				}
			}
			else{  // base == 15

				for (i = 0, j = 0; i < N_base && j < N_base; i++)
				{
					if (j == i)  //just copy, no bit reverse
					{
						unsigned int utemp;
						utemp = i*base;
						tempbuf[utemp] = x[i];
						tempbuf[utemp + 1] = x[N_base + i];
						tempbuf[utemp + 2] = x[N_base * 2 + i];
						tempbuf[utemp + 3] = x[N_base * 3 + i];
						tempbuf[utemp + 4] = x[N_base * 4 + i];
						tempbuf[utemp + 5] = x[N_base * 5 + i];
						tempbuf[utemp + 6] = x[N_base * 6 + i];
						tempbuf[utemp + 7] = x[N_base * 7 + i];
						tempbuf[utemp + 8] = x[N_base * 8 + i];
						tempbuf[utemp + 9] = x[N_base * 9 + i];
						tempbuf[utemp + 10] = x[N_base * 10 + i];
						tempbuf[utemp + 11] = x[N_base * 11 + i];
						tempbuf[utemp + 12] = x[N_base * 12 + i];
						tempbuf[utemp + 13] = x[N_base * 13 + i];
						tempbuf[utemp + 14] = x[N_base * 14 + i];
					}
					if (j > i)  // bit reverse and copy
					{
						unsigned int utemp;
						utemp = j*base;
						tempbuf[utemp] = x[i];
						tempbuf[utemp + 1] = x[N_base + i];
						tempbuf[utemp + 2] = x[N_base * 2 + i];
						tempbuf[utemp + 3] = x[N_base * 3 + i];
						tempbuf[utemp + 4] = x[N_base * 4 + i];
						tempbuf[utemp + 5] = x[N_base * 5 + i];
						tempbuf[utemp + 6] = x[N_base * 6 + i];
						tempbuf[utemp + 7] = x[N_base * 7 + i];
						tempbuf[utemp + 8] = x[N_base * 8 + i];
						tempbuf[utemp + 9] = x[N_base * 9 + i];
						tempbuf[utemp + 10] = x[N_base * 10 + i];
						tempbuf[utemp + 11] = x[N_base * 11 + i];
						tempbuf[utemp + 12] = x[N_base * 12 + i];
						tempbuf[utemp + 13] = x[N_base * 13 + i];
						tempbuf[utemp + 14] = x[N_base * 14 + i];

						utemp = i*base;
						tempbuf[utemp] = x[j];
						tempbuf[utemp + 1] = x[N_base + j];
						tempbuf[utemp + 2] = x[N_base * 2 + j];
						tempbuf[utemp + 3] = x[N_base * 3 + j];
						tempbuf[utemp + 4] = x[N_base * 4 + j];
						tempbuf[utemp + 5] = x[N_base * 5 + j];
						tempbuf[utemp + 6] = x[N_base * 6 + j];
						tempbuf[utemp + 7] = x[N_base * 7 + j];
						tempbuf[utemp + 8] = x[N_base * 8 + j];
						tempbuf[utemp + 9] = x[N_base * 9 + j];
						tempbuf[utemp + 10] = x[N_base * 10 + j];
						tempbuf[utemp + 11] = x[N_base * 11 + j];
						tempbuf[utemp + 12] = x[N_base * 12 + j];
						tempbuf[utemp + 13] = x[N_base * 13 + j];
						tempbuf[utemp + 14] = x[N_base * 14 + j];
					}
					k = N_base / 2;
					while (k >= 2 && j >= k)
					{
						j -= k;
						k >>= 1;
					}
					j += k;
				}
			}
		}

		// perform base transform 4 or 5
		if (base == 4)
		{

			/*************************************************************************
			Forward FFT transform (in-place) for size of 4.

			In order to make the comments clear, "y" and "Y" are used to denote input
			and output sequences. "yr" or "Yr" indicates real part of "y" or "Y", and
			"yi" or "Yi" indicates imaginary part. "x" is used to denote memeory positions.

			For example, input is y[0], y[1], y[2], y[3], loaded from x[i], x[i+2], x[i+1],
			x[i+3]. Output is Y[0], Yr[1], Yr[2], Yi[1], saved into x[i], x[i+1], x[i+2],
			x[i+3], respectively; Note that input is bit-reverse order

			Y[0]  = y[0] + y[1] + y[2] + y[3];
			Yr[1] = y[0] - y[2];
			Yr[2] = y[0] - y[1] + y[2] - y[3];
			Yi[1] = -y[1] + y[3];

			*************************************************************************/
			for (i = 0; i < FFTSize; i += 4)
			{
				float temp1, temp2;
				temp1 = x[i] + x[i + 2] + x[i + 1] + x[i + 3];
				temp2 = x[i] - x[i + 2] + x[i + 1] - x[i + 3];
				x[i + 1] = x[i] - x[i + 1];
				x[i + 3] = x[i + 3] - x[i + 2];
				x[i] = temp1;
				x[i + 2] = temp2;
			}

		}
		else{
			if (base == 5)
			{
				/*********************************************************************************
				forward base transform for size of 5
				input is real, loaded from x2; output is complex, saved into x
				input: xin (real):                 output: xout (complex)
				x[0]                                X[0]
				x[1]                                Xr[1]
				x[2]                                Xr[2]
				x[3]                                Xi[2]
				x[4]                                Xi[1]
				*********************************************************************************/
				for (i = 0; i + 4 < FFTSize; i += 5)
				{
					x[i] = tempbuf[i] + tempbuf[i + 1] + tempbuf[i + 2] + tempbuf[i + 3] + tempbuf[i + 4];
					x[i + 1] = tempbuf[i] + tempbuf[i + 1] * wr1 + tempbuf[i + 2] * wr2 + tempbuf[i + 3] * wr3 + tempbuf[i + 4] * wr4;
					x[i + 4] = +tempbuf[i + 1] * wi1 + tempbuf[i + 2] * wi2 + tempbuf[i + 3] * wi3 + tempbuf[i + 4] * wi4;
					x[i + 2] = tempbuf[i] + tempbuf[i + 1] * wr2 + tempbuf[i + 2] * wr4 + tempbuf[i + 3] * wr1 + tempbuf[i + 4] * wr3;
					x[i + 3] = +tempbuf[i + 1] * wi2 + tempbuf[i + 2] * wi4 + tempbuf[i + 3] * wi1 + tempbuf[i + 4] * wi3;
				}

			}
			else{
				// base transform of size 15
				for (i = 0; i + 14 < FFTSize; i += 15)
					FwdFFT_base15(tempbuf + i, x + i);
			}
		}

		// now we can do in-place butterfly calculation
		step = base;
		while (step<FFTSize)
		{
			// INT note:  the following nested loops can be straighten to a single loop, with
			// if statement in the loop. Need to test which way is faster.
			for (i = 0; i<FFTSize; i += (step * 2))
			{
				float temp;

				// first butterfly
				temp = x[i] - x[i + step];
				x[i] += x[i + step];
				x[i + step] = temp;

				// last butterfly
				// update last butterfly only when step is a even number. Inputs are real, 
				// outputs are complex. Y[p], and Y[q] are complex conjugate to each other,
				// only Y[p] is saved. Yr[p] saved into x[i+step/2], Yi[p] saved into
				// x[i+step+step/2];  
				// Note Yr[p] = y[p] in same place. No need to update.
				//      Yi[p] = -y[q] in same place, only need to take negation.
				if (step % 2 == 0)
					x[i + step + step / 2] = -x[i + step + step / 2];

				// other butterfly. both inputs and outputs are complex
				for (j = 1; j<(step + 1) / 2; j++)
				{
					float wr, wi, tempr, tempi;
					int r, pr, pi, qr, qi;

					pr = i + j;
					pi = i + step - j;
					qr = pr + step;
					qi = pi + step;
					r = j*N_base / 2;      // power of W  (0 <= r <= FFTSize/4)

					wr = cos_tab[-r];
					wi = -sin_tab[r];

					/**************************************************************************
					butterfly computation. this is more complicated than complex butterfly
					* - real values,  o - complex values. the lines indicate values that are
					complex conjugates of each other

					inputs             outputs                   Index    memory
					*                   *                               x[i]
					yr[p] o ----|  --\    /-- o --------| Yr[p]       pr      x[i+1]
					o --| |     \  /    o ------| |                     x[i+2]
					o --| |      \/     o ----| | |                     x[i+3]
					yi[p] o ----|      /\     o --| | | | Yr[q]       pi      x[i+4]
					*           /  \    *   | | | |                     x[i+5]
					yr[q] o ----|  --/    \-- o --| | | | -Yi[q]      qr      x[i+6]
					o --| |             o ----| | |                     x[i+7]
					o --| |             o ------| |                     x[i+8]
					yi[q] o ----|             o --------| Yi[p]       qi      x[i+9]

					Yr[p] = yr[p] + (yr[q]*wr - yi[q]*wi);
					Yi[p] = yi[p] + (yi[q]*wr + yr[q]*wi);
					Yr[q] = yr[p] - (yr[q]*wr - yi[q]*wi);
					Yi[q] = yi[p] - (yi[q]*wr + yr[q]*wi);

					In the following code, tempr = yr[q]*wr - yi[q]*wi
					tempi = yi[q]*wr + yr[q]*wi

					After calculation,  Yr[p] saved to x[pr]
					Yi[p] saved to x[qi]
					Yr[q] saved to x[pi]
					-Yi[q] saved to x[qr]

					**************************************************************************/

					tempr = x[qr] * wr - x[qi] * wi;
					tempi = x[qr] * wi + x[qi] * wr;

					x[qi] = x[pi] + tempi;   //update Yi[p] saved to x[qi]
					x[qr] = -x[pi] + tempi;   //update Yi[q] saved to x[qr]

					x[pi] = x[pr] - tempr;    //update Yr[q] saved to x[pi]
					x[pr] += tempr;           //update Yr[p] saved to x[pr]
				} // inner loop

			} // outer loop

			step *= 2;
			N_base /= 2;     // N_base*step = FFT_SIZE

		} // while loop
		delete[] sin_tab;
		delete[] tempbuf;
	}

	/***************************************************************************
	Function name: AecCcsInvFFT  - CCS Inverse FFT.

	Decription:
	This function supports two output coefficient orders - DFT_COEFF_ORDER_AEC and DFT_COEFF_ORDER_NRM.
	For DFT_COEFF_ORDER_AEC, pfltInData is FFTSize long, otherwise it is FFTSize+2 long.

	pFltSinTable points to a FFTSize/4+1 long sin table which is necessary for ANSI C FFT function AecInvFFT.
	IntelIPP function does not use it. Application is responsible for providing corrent sintable to this function.
	Please refer to filtbank.c for example.

	If FFTSize if power of 2 and DFT_COEFF_ORDER_AEC is specified, this function use in-place computation;
	otherwise a temp buffer will be used internally, but you can still use same buffer for input and
	output sequences. If input and output buffers are different, the input buffer data will not be changed.

	Input:  pXformInfo:  AEC transform structure information
	x:  Pointer of complex frequency sequence, only  0 - FFTSize/2 complex samples
	are saved as:
	X[0], Xr[1], Xr[2], ..., Xr[FFTSize/2-1], X[FFTSize/2], Xi[FFTSize/2-1], ..., Xi[2], Xi[1];
	sin_table:  sin table containing values of sin(2*pi*(0:FFTSize/4)/FFTSize). The size
	of table is FFTSize/4+1;
	FFTSize:  length of time sequence, FFTSize must 2^n, 5*2^n, or 15*2^n
	fCoeffOrder: DFT_COEFF_ORDER_AEC or DFT_COEFF_ORDER_NRM

	Output: x:  Pointer of real time sequence in normal order x[0], x[1], ... , x[FFTSize-1];
	The orignal input is replaced by output.

	Implementation: Based on decimation-in-frequency radix-2 algorithm. The butterfly
	is progressively calculated untill the FFT size reaches 4 or 5.
	Then base transform for size 4 or 5 is calculated directly.

	Note:  Inverse transform is not scaled by 1/FFTSize.

	Reference:  Sorensen, Jones, at el, "Real Valued Fast Fourier Transform
	Algorithms," IEEE Trans. Acoust., Speech, Signal Processing,
	Vol ASSP-36, 1987.

	Qin Li, Feb 18, 2005
	***************************************************************************/
	void FFT::AecCcsInvFFT(float * xin, float * xout, unsigned int FFTSize) {

		// sin and cos table for base 5 FFT
		const float wr1 = 0.309016994374947f;  // cos(2pi/5)
		const float wr2 = -0.809016994374947f;  // cos(4pi/5)
		const float wr3 = -0.809016994374947f;  // cos(6pi/5)
		const float wr4 = 0.309016994374947f;  // cos(8pi/5)
		const float wi1 = 0.951056516295154f;  // sin(2pi/5) 
		const float wi2 = 0.587785252292473f;  // sin(4pi/5)
		const float wi3 = -0.587785252292473f;  // sin(6pi/5)
		const float wi4 = -0.951056516295154f;  // sin(8pi/5)
		float * tempbuf = new float[2 * FFTSize];
		float * sin_tab = new float[FFTSize / 4 + 1];
		for (unsigned int i = 0; i <= FFTSize / 4; i++) {
			sin_tab[i] = (float)sinf(2.0f * (float)PI * i / FFTSize);
		}
		unsigned int base, step, N_base;
		unsigned int i, j, k;   // loop indices
		float * cos_tab, *x;
		cos_tab = sin_tab + FFTSize / 4;
		if (xin != xout)
		{
			memcpy(xout, xin, FFTSize*sizeof(float));
		}
		x = xout;

		if ((FFTSize & (-(int)FFTSize)) == FFTSize) // detect if FFTSize is power of 2
		{
			base = 4;
		}
		else{
			if (FFTSize % 15)
				base = 5;
			else
				base = 15;
		}
		N_base = FFTSize / base;  //number of base transforms

		// Perform butterfly calculation until reach base transform 4 or 5;
		step = FFTSize;
		N_base = 1;
		while (step > base)
		{
			step /= 2;
			N_base *= 2;    // N_base*step = FFT_SIZE
			for (i = 0; i<FFTSize; i += step * 2)
			{
				float temp;

				// first butterfly
				temp = x[i] - x[i + step];
				x[i] += x[i + step];
				x[i + step] = temp;

				// last butterfly
				// update last butterfly only when step is a even number. Inputs are conjugate 
				// symmetric; outputs are real. 
				if (step % 2 == 0)
				{
					x[i + step / 2] *= 2;
					x[i + step + step / 2] = -x[i + step + step / 2] * 2;
				}

				// other butterfly. both inputs and outputs are complex
				for (j = 1; j<(step + 1) / 2; j++)
				{
					int r, pr, pi, qr, qi;
					float wr, wi, tempr, tempi;

					pr = i + j;
					pi = i + step - j;
					qr = pr + step;
					qi = pi + step;
					r = j*N_base / 2;      // power of W  (0 <= r <= FFTSize/4)

					wr = cos_tab[-r];
					wi = sin_tab[r];   // this if for inverse tansform. Take negation for 
					// forward transform.

					/**************************************************************************
					butterfly computation. this is more complicated than complex butterfly
					* - real values,  o - complex values. the lines indicate values that are
					complex conjugates of each other.

					outputs                  inputs              index   memory
					*                       *                          x[i]
					Yr[p]  o --------|  --\    /-- o ----| yr[p]       pr     x[i+1]
					o ------| |     \  /    o --| |                    x[i+2]
					o ----| | |      \/     o --| |                    x[i+3]
					Yr[q]  o --| | | |      /\     o ----| yi[p]       pi     x[i+4]
					*   | | | |     /  \    *                          x[i+5]
					-Yi[q]  o --| | | |  --/    \-- o ----| yr[q]       qr     x[i+6]
					o ----| | |             o --| |                    x[i+7]
					o ------| |             o --| |                    x[i+8]
					Yi[p]  o --------|             o ----| yi[q]       qi     x[i+9]

					(Decimation in frequency)
					y[p] = Y[p] + Y[q];
					y[q] = (Y[p] - Y[q])*W;
					yr[q] = (Yr[p] - Yr[q])*wr - (Yi[p] - Yi[q])*wi;
					yi[q] = (Yr[p] - Yr[q])*wi + (Yi[p] - Yi[q])*wr;

					After calculation,  Yr[p] form x[pr], yr[p] saved to x[pr]
					Yi[p] form x[qi], yi[p] saved to x[pi]
					Yr[q] form x[pi], yr[q] saved to x[qr]
					-Yi[q] form x[qr], yi[q] saved to x[qi]

					**************************************************************************/
					tempr = x[pr] - x[pi];  // tempr = Yr[p] - Yr[q]
					tempi = x[qr] + x[qi];  // tempi = Yi[p] - Yi[q]

					x[pr] += x[pi];         // update yr[p] saved to x[pr]
					x[pi] = x[qi] - x[qr];  // update yi[p] saved to x[pi]

					x[qr] = tempr*wr - tempi*wi; // update yr[q] saved to x[qr]
					x[qi] = tempr*wi + tempi*wr; // update yi[q] saved to x[qi]
				}  // inner loop
			}  // outer loop
		} // while loop

		// Base transform
		if (base == 4)
		{
			/*************************************************************************
			Inverse FFT transform (in-place) for size of 4.

			In order to make the comments clear, "y" and "Y" are used to denote input
			and output sequences. "yr" or "Yr" indicates real part of "y" or "Y", and
			"yi" or "Yi" indicates imaginary part. "x" is used to denote memeory positions.

			For example, iutput is Y[0], Yr[1], Yr[2], Yi[1], loaded from x[i], x[i+1],
			x[i+2], x[i+3]; output is y[0], y[1], y[2], y[3], saved into x[i], x[i+2],
			x[i+1], x[i+3], respectively; Note that ouput is in bit-reverse order.

			y[0] = Y0 + Yr[1] + Yr[2] + Yr[3] = Y0 + Yr[2] + 2*Yr[1];
			y[1] = Y0 - Yi[1] - Yr[2] + Yi[3] = Y0 - Yr[2] - 2*Yi[1];
			y[2] = Y0 - Yr[1] + Yr[2] - Yr[3] = Y0 + Yr[2] - 2*Yr[1];
			y[3] = Y0 + Yi[1] - Yr[2] - Yi[3] = Y0 - Yr[2] + 2*Yi[1];

			*************************************************************************/
			for (i = 0; i < FFTSize; i += 4)
			{
				float temp1, temp2;
				temp1 = x[i] + x[i + 2];     // Y[0] + Yr[2] real
				temp2 = x[i] - x[i + 2];     // Y[0] - Yr[2] real

				x[i] = temp1 + x[i + 1] * 2;  // update y[0], saved into x[i]
				x[i + 2] = temp2 - x[i + 3] * 2;  // update y[1], saved into x[i+2]
				x[i + 1] = temp1 - x[i + 1] * 2;  // update y[2], saved into x[i+1]
				x[i + 3] = temp2 + x[i + 3] * 2;  // update y[3], saved into x[i+3]
			}
		}
		else{
			if (base == 5)
			{
				/********************************************************************************
				forward base transform for size of 5
				input is complex from x, output is real saved into x2. Use seperate input and output
				buffers because cann't do in-place bit reversal for base 5.

				input: xin (conjugate symmetric)       output: xout (real)
				Yr[1]                                y[1]
				Yr[2]                                y[2]
				Yi[2]                                y[3]
				Yi[1]                                y[4]
				********************************************************************************/
				for (i = 0; i + 4 < FFTSize; i += 5)
				{
					tempbuf[i] = x[i] + x[i + 1] * 2 + x[i + 2] * 2;
					// y1 = Y0 + (Yr1*Wr1 - Yi1*Wi1)*2 + (Yr2*Wr2 - Yi2*Wi2)*2;
					tempbuf[i + 1] = x[i] + (x[i + 1] * wr1 - x[i + 4] * wi1 + x[i + 2] * wr2 - x[i + 3] * wi2) * 2;
					// y2 = Y0 + (Yr1*Wr2 - Yi1*Wi2)*2 + (Yr2*Wr4 - Yi2*Wi4)*2;
					tempbuf[i + 2] = x[i] + (x[i + 1] * wr2 - x[i + 4] * wi2 + x[i + 2] * wr4 - x[i + 3] * wi4) * 2;
					// y3 = Y0 + (Yr1*Wr3 - Yi1*Wi3)*2 + (Yr2*Wr1 - Yi2*Wi1)*2;
					tempbuf[i + 3] = x[i] + (x[i + 1] * wr3 - x[i + 4] * wi3 + x[i + 2] * wr1 - x[i + 3] * wi1) * 2;
					// y4 = Y0 + (Yr1*Wr4 - Yi1*Wi4)*2 + (Yr2*Wr3 - Yi2*Wi3)*2;
					tempbuf[i + 4] = x[i] + (x[i + 1] * wr4 - x[i + 4] * wi4 + x[i + 2] * wr3 - x[i + 3] * wi3) * 2;
				}
			}
			else{
				for (i = 0; i + 14 < FFTSize; i += 15)
					InvFFT_base15(x + i, tempbuf + i);
			}
		}

		// convert bit-reversed order to normal order
		if (base == 4)  // FFTSize is power of 2
		{
			for (i = 0, j = 0; i < FFTSize; i++)
			{
				if (j>i)
				{
					float temp;
					temp = x[j];
					x[j] = x[i];
					x[i] = temp;
				}
				k = FFTSize / 2;
				while (k >= 2 && j >= k)
				{
					j -= k;
					k >>= 1;
				}
				j += k;
			}
		}
		else{
			if (base == 5)
			{
				/***********************************************************************
				This is just the inverse process of the bit reversal in the forward fft function.

				For example, for FFT of size 20:
				bit-reversed order              normal order
				0   (00000)                   0   (00000)
				4   (00100)                   1   (00001)
				8   (01000)                   2   (00010)
				12  (01100)                   3   (00011)
				16  (10000)                   4   (00100)
				2   (00010)                   5   (00101)
				6   (00110)                   6   (00110)
				10  (01010)                   7   (00111)
				14  (01110)                   8   (01000)
				18  (10010)                   9   (01001)
				1   (00001)                   10  (01010)
				5   (00101)                   11  (01011)
				9   (01001)                   12  (01100)
				13  (01101)                   13  (01101)
				17  (10001)                   14  (01110)
				3   (00011)                   15  (01111)
				7   (00111)                   16  (10000)
				11  (01011)                   17  (10001)
				15  (01111)                   18  (10010)
				19  (10011)                   19  (10011)
				***********************************************************************/

				for (i = 0, j = 0; i < N_base && j < N_base; i++)
				{
					if (j == i)  //just copy, no bit reverse
					{
						unsigned int utemp;
						utemp = i*base;
						x[i] = tempbuf[utemp];
						x[N_base + i] = tempbuf[utemp + 1];
						x[N_base * 2 + i] = tempbuf[utemp + 2];
						x[N_base * 3 + i] = tempbuf[utemp + 3];
						x[N_base * 4 + i] = tempbuf[utemp + 4];
					}
					if (j > i)  // bit reverse and copy
					{
						unsigned int utemp;
						utemp = j*base;
						x[i] = tempbuf[utemp];
						x[N_base + i] = tempbuf[utemp + 1];
						x[N_base * 2 + i] = tempbuf[utemp + 2];
						x[N_base * 3 + i] = tempbuf[utemp + 3];
						x[N_base * 4 + i] = tempbuf[utemp + 4];

						utemp = i*base;
						x[j] = tempbuf[utemp];
						x[N_base + j] = tempbuf[utemp + 1];
						x[N_base * 2 + j] = tempbuf[utemp + 2];
						x[N_base * 3 + j] = tempbuf[utemp + 3];
						x[N_base * 4 + j] = tempbuf[utemp + 4];
					}
					k = N_base / 2;
					while (k >= 2 && j >= k)
					{
						j -= k;
						k >>= 1;
					}
					j += k;
				}
			}
			else{  // base = 15

				for (i = 0, j = 0; i < N_base && j < N_base; i++)
				{
					if (j == i)  //just copy, no bit reverse
					{
						unsigned int utemp;
						utemp = i*base;
						x[i] = tempbuf[utemp];
						x[N_base + i] = tempbuf[utemp + 1];
						x[N_base * 2 + i] = tempbuf[utemp + 2];
						x[N_base * 3 + i] = tempbuf[utemp + 3];
						x[N_base * 4 + i] = tempbuf[utemp + 4];
						x[N_base * 5 + i] = tempbuf[utemp + 5];
						x[N_base * 6 + i] = tempbuf[utemp + 6];
						x[N_base * 7 + i] = tempbuf[utemp + 7];
						x[N_base * 8 + i] = tempbuf[utemp + 8];
						x[N_base * 9 + i] = tempbuf[utemp + 9];
						x[N_base * 10 + i] = tempbuf[utemp + 10];
						x[N_base * 11 + i] = tempbuf[utemp + 11];
						x[N_base * 12 + i] = tempbuf[utemp + 12];
						x[N_base * 13 + i] = tempbuf[utemp + 13];
						x[N_base * 14 + i] = tempbuf[utemp + 14];
					}
					if (j > i)  // bit reverse and copy
					{
						unsigned int utemp;
						utemp = j*base;
						x[i] = tempbuf[utemp];
						x[N_base + i] = tempbuf[utemp + 1];
						x[N_base * 2 + i] = tempbuf[utemp + 2];
						x[N_base * 3 + i] = tempbuf[utemp + 3];
						x[N_base * 4 + i] = tempbuf[utemp + 4];
						x[N_base * 5 + i] = tempbuf[utemp + 5];
						x[N_base * 6 + i] = tempbuf[utemp + 6];
						x[N_base * 7 + i] = tempbuf[utemp + 7];
						x[N_base * 8 + i] = tempbuf[utemp + 8];
						x[N_base * 9 + i] = tempbuf[utemp + 9];
						x[N_base * 10 + i] = tempbuf[utemp + 10];
						x[N_base * 11 + i] = tempbuf[utemp + 11];
						x[N_base * 12 + i] = tempbuf[utemp + 12];
						x[N_base * 13 + i] = tempbuf[utemp + 13];
						x[N_base * 14 + i] = tempbuf[utemp + 14];

						utemp = i*base;
						x[j] = tempbuf[utemp];
						x[N_base + j] = tempbuf[utemp + 1];
						x[N_base * 2 + j] = tempbuf[utemp + 2];
						x[N_base * 3 + j] = tempbuf[utemp + 3];
						x[N_base * 4 + j] = tempbuf[utemp + 4];
						x[N_base * 5 + j] = tempbuf[utemp + 5];
						x[N_base * 6 + j] = tempbuf[utemp + 6];
						x[N_base * 7 + j] = tempbuf[utemp + 7];
						x[N_base * 8 + j] = tempbuf[utemp + 8];
						x[N_base * 9 + j] = tempbuf[utemp + 9];
						x[N_base * 10 + j] = tempbuf[utemp + 10];
						x[N_base * 11 + j] = tempbuf[utemp + 11];
						x[N_base * 12 + j] = tempbuf[utemp + 12];
						x[N_base * 13 + j] = tempbuf[utemp + 13];
						x[N_base * 14 + j] = tempbuf[utemp + 14];
					}
					k = N_base / 2;
					while (k >= 2 && j >= k)
					{
						j -= k;
						k >>= 1;
					}
					j += k;
				}
			}
		}
	}
}

