#ifndef MCLT_H_
#define MCLT_H_

#include <complex>
#include <memory>
#include "Dsplib.h"
#include "fftw3.h"
#include "Math.h"

class MCLT {
public:
	/// constructor with the size.
	MCLT(int M);
	~MCLT();
	/// compute MCLT.
	/// the input size is 2 * M and the output are M complex numbers.
	void compute(double* input, double* re, double* im);
private:
	int m_M;
	double m_coeff;
	std::unique_ptr<double[]> m_ha;
	std::unique_ptr<double[]> m_u;
	std::unique_ptr<double[]> m_v;
};

#endif /* MCLT_H_ */
