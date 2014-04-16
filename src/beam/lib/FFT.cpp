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
		for (int i = 0; i < FRAME_SIZE; ++i){
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
		fftwf_plan p = fftwf_plan_dft_c2r_1d(2 * FRAME_SIZE, in, m_current, FFTW_ESTIMATE);
		fftwf_execute(p);
		fftwf_destroy_plan(p);

		int two_frame_size = 2 * FRAME_SIZE;
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
			memcpy_s(xout, FFTSize*sizeof(float), xin, FFTSize*sizeof(float));
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
}

