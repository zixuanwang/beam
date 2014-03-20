#ifndef NOISESUPPRESSOR_H_
#define NOISESUPPRESSOR_H_

#include "Utils.h"

namespace Beam{
	class NoiseSuppressor {
	public:
		NoiseSuppressor();
		NoiseSuppressor(float frequency, int frame_size, float adaptive_tau, float suppress);
		~NoiseSuppressor();
		void init(float frequency, int frame_size, float adaptive_tau, float suppress);
		void phase_compensation(std::vector<std::complex<float> >& output);
		void noise_compensation(std::vector<std::complex<float> >& output);
		void set_suppress(float suppress);
	private:
		float m_frame_duration;
		float m_phase_adaptive_tau;
		float m_speed_adaptive_tau;
		float m_noise_adaptive_tau;
		float m_sampling_rate;
		// phase compensation data
		int m_phase_num_frames;
		std::vector<std::complex<float> > m_phase_model;
		std::vector<float> m_phase_model_variance;
		std::vector<std::complex<float> > m_phase_model_update;
		std::vector<std::complex<float> > m_phase_model_prev;
		// noise compensation data
		int m_noise_num_frames;
		std::vector<float> m_noise_model;
		std::vector<float> m_noise_model_variance;
		float m_suppress;
	};
}

#endif /* NOISESUPPRESSOR_H_ */
