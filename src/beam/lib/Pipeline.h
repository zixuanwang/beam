#ifndef PIPELINE_H_
#define PIPELINE_H_

#include <iostream>
#include "Beamformer.h"
#include "Calibrator.h"
#include "DelaySumBeamformer.h"
#include "DeReverb.h"
#include "DSPFilter.h"
#include "FFT.h"
#include "GlobalConfig.h"
#include "GSCBeamformer.h"
#include "MCLT.h"
#include "MsrNS.h"
#include "NoiseSuppressor.h"
#include "SoundSourceLocalizer.h"
#include "Tracker.h"
#include "WavReader.h"
#include "WavWriter.h"

namespace Beam{
	class Pipeline{
	public:
		static Pipeline* instance();
		void phase_compensation(float* fft_ptr, bool analysis);
		/// input should have FRAME_SIZE.
		void convert_input(std::vector<std::complex<float> >& input, float* fft_ptr);
		/// fft_ptr should have 2 * FRAME_SIZE.
		void convert_output(const std::vector<std::complex<float> >& output, float* fft_ptr);
		/// process frame.
		void process(float input[MAX_MICROPHONES][FRAME_SIZE], float output[FRAME_SIZE]);
		/// noise suppression.
		void suppress_noise(float* fft_ptr);

		// multi-channel inputs.
		void preprocess(std::vector<std::complex<float> >* input);
		// if the sound source can be localized, the angle is store in p_angle.
		void source_localize(std::vector<std::complex<float> >* input, float* p_angle);
		void dereverbration(std::vector<std::complex<float> >* input);
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
		NoiseSuppressor m_pre_suppressor[MAX_MICROPHONES]; // for noise suppression.
		DeReverb m_dereverb[MAX_MICROPHONES];
		NoiseSuppressor m_out_noise_suppressor; // for frequency shifting the output.
		Tracker m_noise_floor; // VAD
		SoundSourceLocalizer m_ssl; // SSL
		Calibrator m_calibrator; // Calibrator
		Beamformer m_beamformer; // BF
		float m_confidence;
		float m_angle; // sound source angle
		bool m_voice_found; // result of VAD
		bool m_source_found;
		int m_frame_number;
		// gains.
		std::vector<std::complex<float> > m_dynamic_gains[MAX_MICROPHONES];
		std::complex<float> m_persistent_gains[MAX_MICROPHONES][MAX_GAIN_SUBBANDS];
		std::vector<std::complex<float> > m_input_channels[MAX_MICROPHONES];
		int m_refresh_gain;
		float m_input_prev[MAX_MICROPHONES][FRAME_SIZE];
		float m_input[MAX_MICROPHONES][2 * FRAME_SIZE];
		float m_output_prev[FRAME_SIZE];
		float m_output[2 * FRAME_SIZE];
		std::vector<std::complex<float> > m_frequency_input[MAX_MICROPHONES];
		std::vector<std::complex<float> > m_frequency_output;
		// timer. every time preprocess is called, the time is updated.
		double m_time;
		MsrNS m_ns;
		MsrVAD m_vad;
	};
}

#endif /* PIPELINE_H_ */
