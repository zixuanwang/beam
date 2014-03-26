#ifndef PIPELINE_H_
#define PIPELINE_H_

#include <iostream>
#include "DSPFilter.h"
#include "MicArrayDescriptor.h"
#include "MicArrayWeights.h"
#include "NoiseSuppressor.h"
#include "SoundSourceLocalizer.h"
#include "Tracker.h"

namespace Beam{
	class Pipeline{
	public:
		static Pipeline* instance();
		void load_profile();
		// multi-channel inputs.
		void preprocess(std::vector<std::complex<float> >* input);
		// if the sound source can be localized, true is returned and angle is store in p_angle.
		// otherwise, false is returned.
		bool source_localize(std::vector<std::complex<float> >* input, double time, float* p_angle);
	private:
		// singleton.
		Pipeline();
		Pipeline(Pipeline&);
		Pipeline& operator=(Pipeline&);
		static Pipeline* p_instance;
		NoiseSuppressor m_pre_noise_suppressor[MAX_MICROPHONES]; // for phase compensation.
		NoiseSuppressor m_ssl_noise_suppressor[MAX_MICROPHONES]; // for noise suppression in ssl.
		std::vector<std::complex<float> > m_ssl_band_pass_filter;
		Tracker m_noise_floor; // VAD
		SoundSourceLocalizer m_ssl;
	};
}

#endif /* PIPELINE_H_ */
