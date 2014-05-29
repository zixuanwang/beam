#include "DeReverb.h"

namespace Beam{
	DeReverb::DeReverb() : m_voice_frame_count(0), m_voice_found(false), m_tail_found(false), m_tail_count(0){
		m_cepstral_mean.assign(FRAME_SIZE, std::complex<float>(0.f, 0.f));
		m_init_energy.assign(FRAME_SIZE, 0.f);
		for (int i = 0; i < FRAME_SIZE; ++i){
			m_energy[i].assign(TAIL_FRAME_SIZE, 0.f);
		}
		m_tau.assign(FRAME_SIZE, 0.02f);
		m_energy_list.assign(FRAME_SIZE, std::list<float>());
	}

	DeReverb::~DeReverb(){
	
	}

	void DeReverb::normalize_cepstral(std::vector<std::complex<float> >& input, bool voice_found){
		if (voice_found){
			for (int bin = 0; bin < FRAME_SIZE; ++bin){
				std::complex<float> cepstral;
				float scale = Utils::abs_complex(input[bin]);
				if (scale > 0.f){
					cepstral.real(logf(scale));
					cepstral.imag(std::arg(input[bin]));
					if (m_voice_frame_count == 0){
						m_cepstral_mean[bin] = cepstral;
					}
					else{
						m_cepstral_mean[bin] = m_cepstral_mean[bin] * 0.95f + cepstral * 0.05f; // running average
					}
					std::complex<float> normalized_cepstral = cepstral - m_cepstral_mean[bin];
					float scale = expf(normalized_cepstral.real());
					float angle = normalized_cepstral.imag();
					input[bin].real(scale * cosf(angle));
					input[bin].imag(scale * sinf(angle));
				}
			}
			++m_voice_frame_count;
		}
	}

	void DeReverb::suppress(std::vector<std::complex<float> >& input, bool voice_found){
		/*if (voice_found){
			for (int bin = 0; bin < FRAME_SIZE; ++bin){
				m_init_energy[bin] = Utils::norm_complex(input[bin]);
			}
			m_tail_found = false;
		}
		if (m_voice_found && !voice_found){
			m_tail_found = true;
			m_tail_count = 0;
		}
		if (m_tail_found){
			for (int bin = 0; bin < FRAME_SIZE; ++bin){
				m_energy[bin][m_tail_count] = Utils::norm_complex(input[bin]);
			}
			++m_tail_count;
			if (m_tail_count >= TAIL_FRAME_SIZE){
				std::vector<float> tau(FRAME_SIZE, 0.f);
				for (int bin = 0; bin < FRAME_SIZE; ++bin){
					tau[bin] = compute_tau(m_init_energy[bin], m_energy[bin]);
				}
				update_tau(tau);
				m_tail_found = false;
			}
		}
		m_voice_found = voice_found;*/
		for (int bin = 0; bin < FRAME_SIZE; ++bin){
			float element = Utils::abs_complex(input[bin]);
			if (m_energy_list[bin].size() == TAIL_FRAME_SIZE){
				float model = 0.f;
				int i = 0;
				float t = (float)FRAME_SIZE / SAMPLE_RATE;
				for (auto energy : m_energy_list[bin]){
					model += energy * expf(-1.f / 0.02f * (TAIL_FRAME_SIZE - i) * t);
					++i;
				}
				float gain = 0.f;
				if (element > model && voice_found){
					gain = (element * element - 0.9f * model * model) / (element * element);
				}
				else{
					gain = 1.f - 0.9f;
				}
				input[i] *= gain;
			}
			if (m_energy_list[bin].size() >= TAIL_FRAME_SIZE){
				m_energy_list[bin].pop_front();
			}
			m_energy_list[bin].push_back(element);
		}
	}

	float DeReverb::compute_tau(float init_energy, const std::vector<float>& tail_energy){
		if (init_energy == 0.f){
			return 0.02f;
		}
		float t = (float)FRAME_SIZE / SAMPLE_RATE;
		float tau_inv = 0.f;
		int valid_count = 0;
		float log_init_energy = logf(init_energy);
		for (int i = 0; i < TAIL_FRAME_SIZE; ++i){
			if (tail_energy[i] != 0.f){
				float log_diff = log_init_energy - logf(tail_energy[i]);
				if (log_diff > 0.f){
					tau_inv += (log_init_energy - logf(tail_energy[i])) / (t * (float)(i + 1));
					++valid_count;
				}
			}
		}
		if (tau_inv != 0.f && valid_count > 0){
			return (float)valid_count / tau_inv;
		}
		return 0.02f;
	}

	void DeReverb::update_tau(const std::vector<float>& tau){
		for (int bin = 0; bin < FRAME_SIZE; ++bin){
			if (m_tau[bin] == 0.f){
				m_tau[bin] = tau[bin];
			}
			else{
				m_tau[bin] = 0.8f * m_tau[bin] + 0.2f * tau[bin];
			}
		}
	}
}