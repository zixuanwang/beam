#ifndef MCLT_H_
#define MCLT_H_

#include "fftw3.h"
#include "Utils.h"

namespace Beam{
	class MCLT {
	public:
		/// constructor with the size.
		MCLT(int M);
		~MCLT();
		/// compute MCLT.
		/// the input size is 2 * M and the output are M complex numbers.
		void compute(const std::vector<float>& input, std::vector<std::complex<float> >& output);
	private:
		int m_M;
		double m_coeff;
		std::unique_ptr<double[]> m_ha;
		std::unique_ptr<double[]> m_u;
		std::unique_ptr<double[]> m_v;
	};
}

#endif /* MCLT_H_ */
