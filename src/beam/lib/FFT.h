#ifndef FFT_H_
#define FFT_H_

#include "fftw3.h"
#include "GlobalConfig.h"
#include "Utils.h"

namespace Beam{
	class FFT {
	public:
		/// constructor with the size.
		FFT(int n);
		~FFT();
		/// compute FFT.
		void compute(const std::vector<float>& input, std::vector<std::complex<float> >& output);
	private:
		int m_n;
		double m_ha[FRAME_SIZE];
	};
}

#endif /* FFT_H_ */
