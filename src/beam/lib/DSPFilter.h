#ifndef DSPFILTER_H_
#define DSPFILTER_H_

#include "GlobalConfig.h"
#include "Utils.h"

namespace Beam{
	class DSPFilter
	{
	public:
		DSPFilter();
		~DSPFilter();
		/// filters using MCLT.
		static void low_pass_mclt(std::vector<std::complex<float> >& spectrum, float pass_freq, float stop_freq);
		static void high_pass_mclt(std::vector<std::complex<float> >& spectrum, float pass_freq, float stop_freq);
		static void band_pass_mclt(std::vector<std::complex<float> >& spectrum, float low_stop_freq, float low_pass_freq, float high_pass_freq, float high_stop_freq);
		/// filters using FFT.
		static void low_pass_fft(std::vector<std::complex<float> >& spectrum, float pass_freq, float stop_freq);
		static void high_pass_fft(std::vector<std::complex<float> >& spectrum, float pass_freq, float stop_freq);
		static void band_pass_fft(std::vector<std::complex<float> >& spectrum, float low_stop_freq, float low_pass_freq, float high_pass_freq, float high_stop_freq);

	};
}

#endif /* DSPFILTER_H_ */
