#ifndef DSPFILTER_H_
#define DSPFILTER_H_

#include "GlobalConfig.h"
#include "MCLT.h"

namespace Beam{
	class DSPFilter
	{
	public:
		DSPFilter();
		~DSPFilter();
		static void low_pass(std::vector<std::complex<float> >& spectrum, float pass_freq, float stop_freq);
		static void high_pass(std::vector<std::complex<float> >& spectrum, float pass_freq, float stop_freq);
		static void band_pass(std::vector<std::complex<float> >& spectrum, float low_stop_freq, float low_pass_freq, float high_pass_freq, float high_stop_freq);

	};
}

#endif /* DSPFILTER_H_ */
