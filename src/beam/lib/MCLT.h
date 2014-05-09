#ifndef MCLT_H_
#define MCLT_H_

#include "fftw3.h"
#include "GlobalConfig.h"
#include "Utils.h"

namespace Beam{
	class MCLT {
	public:
		/// constructor with the size.
		MCLT();
		~MCLT();
		/// compute MCLT. make sure that output is allocated.
		/// the input must have 2 * FRAME_SIZE.
		/// the output must have FRAME_SIZE.
		void analyze(std::vector<float>& input, std::vector<std::complex<float> >& output);
		/// compute IMCLT. make sure that output is allocated.
		/// the input must have FRAME_SIZE.
		/// the output must have 2 * FRAME_SIZE.
		/// the output has half frame delay.
		void synthesize(std::vector<std::complex<float> >& input, std::vector<float>& output);
	private:
		float m_coeff;
		float m_input[FRAME_SIZE * 2];
		float m_ha[FRAME_SIZE * 2];
		float m_prev_prev[FRAME_SIZE * 2];
		float m_prev[FRAME_SIZE * 2];
		float m_current[FRAME_SIZE * 2];
		float m_u[FRAME_SIZE];
		float m_v[FRAME_SIZE];
	};
}

#endif /* MCLT_H_ */
