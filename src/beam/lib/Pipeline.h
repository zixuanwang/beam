#ifndef PIPELINE_H_
#define PIPELINE_H_

#include <iostream>
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
		// if the sound source can be localized, true is returned and angle is store in p_angle.
		// otherwise, false is returned.
		bool source_localize(std::vector<std::complex<float> >* input, double time, float* p_angle);
		void beamformer(std::vector<std::complex<float> >* input, std::vector<std::complex<float> >& output, double time);
		float smart_calibration(float sound_source, std::vector<std::complex<float> >* input);
		void expand_gain();
		void ansi_bf_msr_process_quad_loop_fast(std::complex<float>* wo0, std::complex<float>* wo1, std::complex<float>* wo2, std::complex<float>* wo3, std::complex<float>& m0, std::complex<float>& m1, std::complex<float>& m2, std::complex<float>& m3, std::complex<float>& w0, std::complex<float>& w1, std::complex<float>& w2, std::complex<float>& w3, float nu, float mu);
	private:
		// singleton.
		Pipeline();
		Pipeline(Pipeline&);
		Pipeline& operator=(Pipeline&);
		static Pipeline* p_instance;
		NoiseSuppressor m_pre_noise_suppressor[MAX_MICROPHONES]; // for phase compensation in the preprocessing.
		NoiseSuppressor m_ssl_noise_suppressor[MAX_MICROPHONES]; // for noise suppression in ssl.
		NoiseSuppressor m_out_noise_suppressor[MAX_MICROPHONES]; // for frequency shifting the output.
		std::vector<std::complex<float> > m_ssl_band_pass_filter;
		Tracker m_noise_floor; // VAD
		SoundSourceLocalizer m_ssl;
		// for calibration
		float m_coordinates[MAX_MICROPHONES];
		float m_coeff[2];
		std::vector<std::complex<float> > m_frequency_filter[MAX_GAIN_SUBBANDS];
		std::vector<std::complex<float> > m_working_frequency;
		int m_refresh_ini;

		// for beamformer
		float m_confidence;
		float m_source_position;
		int m_beam;
		int m_first_bin;
		int m_last_bin;

		//
		std::vector<std::complex<float> > m_pcm_weights[MAX_BEAMS][MAX_MICROPHONES];
		std::vector<std::complex<float> > m_dynamic_gains[MAX_MICROPHONES];
		std::complex<float> m_persistent_gains[MAX_MICROPHONES][MAX_GAIN_SUBBANDS];
	};
}

#endif /* PIPELINE_H_ */
