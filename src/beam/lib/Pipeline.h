#ifndef PIPELINE_H_
#define PIPELINE_H_

#include "DSPFilter.h"
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
		void source_localize(std::vector<std::complex<float> >* input, std::vector<std::complex<float> >& output, double time);
	private:
		// singleton.
		Pipeline();
		Pipeline(Pipeline&);
		Pipeline& operator=(Pipeline&);
		static Pipeline* p_instance;
		NoiseSuppressor m_pre_noise_suppressor[MAX_MICROPHONES]; // for phase compensation.
		NoiseSuppressor m_ssl_noise_suppressor[MAX_MICROPHONES]; // for noise suppression in ssl.
		std::vector<std::complex<float> > m_ssl_band_pass_filter;
	};
}

#endif /* PIPELINE_H_ */
