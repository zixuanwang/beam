#include "DSPFilter.h"

namespace Beam{
	DSPFilter::DSPFilter()
	{
	}

	DSPFilter::~DSPFilter()
	{
	}

	void DSPFilter::low_pass_mclt(std::vector<std::complex<float> >& spectrum, float pass_freq, float stop_freq){
		int bins = (int)spectrum.size();
		int low_index = (int)(pass_freq * bins * 2.f);
		if (low_index <= 0){
			low_index = 1;
		}
		int up_index = (int)(stop_freq * bins * 2.f + 0.5f);
		if (up_index > bins - 1){
			up_index = bins - 1;
		}
		int index = 0;
		//  Transition band - apply cos() shaped slope
		for (index = low_index; index <= up_index; ++index){
			float arg = (((float)index) / bins / 2.f - pass_freq) / (stop_freq - pass_freq);
			float mult = (cosf((float)(arg * PI)) + 1.f) / 2.f;
			spectrum[index] *= mult;
		}
		//  Suppress band - all zeros
		for (; index < bins; ++index){
			spectrum[index] *= 0.f;
		}
	}

	void DSPFilter::high_pass_mclt(std::vector<std::complex<float> >& spectrum, float pass_freq, float stop_freq){
		int bins = (int)spectrum.size();
		int low_index = (int)(stop_freq * bins * 2.f);
		if (low_index <= 0){
			low_index = 1;
		}
		int up_index = (int)(pass_freq * bins * 2.f + 0.5f);
		if (up_index > bins - 1){
			up_index = bins - 1;
		}
		for (int index = 0; index < low_index; ++index){
			spectrum[index] *= 0.f;
		}
		for (int index = low_index; index <= up_index; ++index){
			float arg = (((float)index) / bins / 2.f - stop_freq) / (pass_freq - stop_freq);
			float mult = (1.f - cosf((float)PI * arg)) / 2.f;
			spectrum[index] *= mult;
		}
	}

	void DSPFilter::band_pass_mclt(std::vector<std::complex<float> >& spectrum, float low_stop_freq, float low_pass_freq, float high_pass_freq, float high_stop_freq){
		spectrum.assign(FRAME_SIZE, std::complex<float>(1.f, 0.f));
		low_pass_mclt(spectrum, high_pass_freq, high_stop_freq);
		high_pass_mclt(spectrum, low_pass_freq, low_stop_freq);
	}

	void DSPFilter::low_pass_fft(std::vector<std::complex<float> >& spectrum, float pass_freq, float stop_freq){
		int bins = (int)spectrum.size();
		int low_index = (int)(pass_freq * bins);
		if (low_index <= 0){
			low_index = 1;
		}
		int up_index = (int)(stop_freq * bins + 0.5f);
		if (up_index > bins / 2){
			up_index = bins / 2;
		}
		//  Transition band - apply cos() shaped slope
		int index = 0;
		for (index = low_index; index <= up_index; ++index){
			float arg = (((float)index) / bins - pass_freq) / (stop_freq - pass_freq);
			float mult = (cosf((float)(arg * PI)) + 1.f) / 2.f;
			spectrum[index] *= mult;
			spectrum[bins - index] *= mult;
		}
		//  Suppress band - all zeros
		for (; index < bins / 2; ++index){
			spectrum[index] *= 0.f;
			spectrum[bins - index] *= 0.f;
		}
	}

	void DSPFilter::high_pass_fft(std::vector<std::complex<float> >& spectrum, float pass_freq, float stop_freq){
		int bins = (int)spectrum.size();
		int low_index = (int)(stop_freq * bins);
		if (low_index <= 0){
			low_index = 1;
		}
		int up_index = (int)(pass_freq * bins + 0.5f);
		if (up_index > bins / 2){
			up_index = bins / 2;
		}
		spectrum[0] *= 0.f;
		for (int index = 1; index < low_index; ++index){
			spectrum[index] *= 0.f;
			spectrum[bins - index] *= 0.f;
		}
		for (int index = low_index; index <= up_index; ++index){
			float arg = (((float)index) / bins - stop_freq) / (pass_freq - stop_freq);
			float mult = (1.f - cosf((float)PI * arg)) / 2.f;
			spectrum[index] *= mult;
			spectrum[bins - index] *= mult;
		}
	}

	void DSPFilter::band_pass_fft(std::vector<std::complex<float> >& spectrum, float low_stop_freq, float low_pass_freq, float high_pass_freq, float high_stop_freq){
		spectrum.assign(FRAME_SIZE, std::complex<float>(1.f, 0.f));
		low_pass_mclt(spectrum, high_pass_freq, high_stop_freq);
		high_pass_mclt(spectrum, low_pass_freq, low_stop_freq);
	}
}
