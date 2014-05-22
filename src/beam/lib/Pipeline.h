#ifndef PIPELINE_H_
#define PIPELINE_H_

#include <iostream>
#include "Beamformer.h"
#include "Calibrator.h"
#include "DelaySumBeamformer.h"
#include "DSPFilter.h"
#include "FFT.h"
#include "GlobalConfig.h"
#include "NoiseSuppressor.h"
#include "SoundSourceLocalizer.h"
#include "Tracker.h"
#include "WavReader.h"
#include "WavWriter.h"

namespace Beam{
	class Pipeline{
	public:
		static Pipeline* instance();
		// multi-channel inputs.
		void preprocess(std::vector<std::complex<float> >* input);
		// if the sound source can be localized, the angle is store in p_angle.
		void source_localize(std::vector<std::complex<float> >* input, float* p_angle);
		void smart_calibration(std::vector<std::complex<float> >* input);
		void beamforming(std::vector<std::complex<float> >* input, std::vector<std::complex<float> >& output);
		void postprocessing(std::vector<std::complex<float> >& input);
		void expand_gain();
	private:
		// singleton.
		Pipeline();
		Pipeline(Pipeline&);
		Pipeline& operator=(Pipeline&);
		static Pipeline* p_instance;
		// components.
		std::vector<std::complex<float> > m_band_pass_filter; // band pass filter
		NoiseSuppressor m_pre_noise_suppressor[MAX_MICROPHONES]; // for phase compensation in the preprocessing.
		NoiseSuppressor m_ssl_noise_suppressor[MAX_MICROPHONES]; // for noise suppression in ssl.
		NoiseSuppressor m_suppressor[MAX_MICROPHONES]; // for testing.
		NoiseSuppressor m_out_noise_suppressor; // for frequency shifting the output.
		Tracker m_noise_floor; // VAD
		SoundSourceLocalizer m_ssl; // SSL
		Calibrator m_calibrator; // Calibrator
		Beamformer m_beamformer; // BF
		float m_confidence;
		float m_angle; // sound source angle
		bool m_source_found;
		// gains.
		std::vector<std::complex<float> > m_dynamic_gains[MAX_MICROPHONES];
		std::complex<float> m_persistent_gains[MAX_MICROPHONES][MAX_GAIN_SUBBANDS];
		std::vector<std::complex<float> > m_input_channels[MAX_MICROPHONES];
		int m_refresh_gain;
		// timer. every time preprocess is called, the time is updated.
		double m_time;
	};
}

#endif /* PIPELINE_H_ */
