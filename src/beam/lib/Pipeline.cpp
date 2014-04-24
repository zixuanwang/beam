#include "Pipeline.h"

namespace Beam{
	Pipeline* Pipeline::p_instance = NULL;

	Pipeline::Pipeline() : m_noise_floor(20.0, 0.04, 30000.0, 0.0){
		// initialize band pass filter.
		DSPFilter::band_pass_mclt(m_band_pass_filter, 500.0 / SAMPLE_RATE, 1000.0 / SAMPLE_RATE, 2000.0 / SAMPLE_RATE, 3500.0 / SAMPLE_RATE);
		// initialize noise suppressors.
		for (int channel = 0; channel < MAX_MICROPHONES; ++channel){
			m_ssl_noise_suppressor[channel].init(SAMPLE_RATE, FRAME_SIZE, 1.f, 10.f);
			m_pre_noise_suppressor[channel].init(SAMPLE_RATE, FRAME_SIZE, 1.f, 1.f);
		}
		m_out_noise_suppressor.init(SAMPLE_RATE, FRAME_SIZE, 1.f, 1.f);
		m_ssl.init(SAMPLE_RATE, FRAME_SIZE);
		// initialize persistent and dynamic gains
		for (int channel = 0; channel < MAX_MICROPHONES; ++channel){
			m_dynamic_gains[channel].assign(FRAME_SIZE, std::complex<float>(1.f, 0.f));
		}
		for (int channel = 0; channel < MAX_MICROPHONES; ++channel){
			for (int sub = 0; sub < MAX_GAIN_SUBBANDS; ++sub){
				m_persistent_gains[channel][sub] = std::complex<float>(1.f, 0.f);
			}
		}
		// initialize m_input_channels
		for (int channel = 0; channel < MAX_MICROPHONES; ++channel){
			m_input_channels[channel].assign(FRAME_SIZE, std::complex<float>(0.f, 0.f));
		}
		// initialize gains.
		expand_gain();
		m_refresh_gain = 0;
		// initialize m_beamformer.
		m_beamformer.init();
		m_confidence = 0.f;
		// initialize m_time.
		m_time = 0.0;
		// initialize m_angle.
		m_angle = 0.f;
		// initialize m_source_found;
		m_source_found = false;
	}

	Pipeline* Pipeline::instance(){
		if (p_instance == NULL){
			p_instance = new Pipeline;
		}
		return p_instance;
	}

	void Pipeline::preprocess(std::vector<std::complex<float> >* input){
		for (int channel = 0; channel < MAX_MICROPHONES; ++channel){
			/*
			for (int bin = 0; bin < FRAME_SIZE; ++bin){
				m_input_channels[channel][bin] *= m_dynamic_gains[channel][bin];
			}
			*/
			m_pre_noise_suppressor[channel].phase_compensation(input[channel]);
		}
		m_time += (double)FRAME_SIZE / (double)SAMPLE_RATE;
	}

	void Pipeline::source_localize(std::vector<std::complex<float> >* input, float* p_angle){
		m_source_found = false;
		//  Apply the SSL band pass filter to the input channels
		//  and have a separate copy of the input channels 
		//  for SSL purposes only
		for (int channel = 0; channel < MAX_MICROPHONES; ++channel){
			for (size_t bin = 0; bin < m_band_pass_filter.size(); ++bin){
				m_input_channels[channel][bin] = input[channel][bin] * m_band_pass_filter[bin];
			}
		}
		//  Noise suppression
		//  We do heavy noise suppression as we don't care about the musical noises
		//  but we do cary to suppress stationaty noises
		double energy = 0.0;
		for (int channel = 0; channel < MAX_MICROPHONES; ++channel){
			m_ssl_noise_suppressor[channel].noise_compensation(m_input_channels[channel]);
			energy += (double)Utils::computeRMS(m_input_channels[channel]);
		}
		energy /= MAX_MICROPHONES;
		double floor = m_noise_floor.nextLevel(m_time, energy);
		if (energy > SSL_RELATIVE_ENERGY_THRESHOLD * floor && energy > SSL_ABSOLUTE_ENERGY_THRESHOLD){
			// sound signal
			float angle;
			float weight;
			m_ssl.process(m_input_channels, &angle, &weight);
			if (weight > SSL_CONTRAST_THRESHOLD){
				m_ssl.process_next_sample(m_time, angle, weight);
				m_source_found = true;
			}
		}
		float std_dev;
		int valid;
		int num;
		m_ssl.get_average(m_time, p_angle, &m_confidence, &std_dev, &num, &valid);
		m_angle = *p_angle;
	}

	void Pipeline::smart_calibration(){
		if (m_source_found){
			float sigma = m_calibrator.calibrate(m_angle, m_input_channels, m_persistent_gains);
			++m_refresh_gain;
			if (m_refresh_gain > 200){
				expand_gain();
				m_refresh_gain = 0;
			}
		}
	}

	void Pipeline::beamforming(std::vector<std::complex<float> >* input, std::vector<std::complex<float> >& output){
		m_beamformer.compute(input, output, m_angle, m_confidence, m_time);
	}

	void Pipeline::postprocessing(std::vector<std::complex<float> >& input){
		m_out_noise_suppressor.frequency_shifting(input);
	}

	void Pipeline::expand_gain(){
		std::complex<float> zero(0.f, 0.f);
		// use MCLT
		float freq_step = (float)SAMPLE_RATE / FRAME_SIZE / 2.f;
		float freq_beg = 0.f;
		if (USE_MCLT){
			freq_beg = freq_step / 2.f;
		}
		for (int index = 0; index < FRAME_SIZE; ++index){
			int interp_high = 1;
			int interp_low = 0;
			float freq = freq_beg + index * freq_step;
			if (freq > KinectConfig::frequency_bands[MAX_GAIN_SUBBANDS - 1][0]){
				freq = KinectConfig::frequency_bands[MAX_GAIN_SUBBANDS - 1][0];
			}
			while (freq >= KinectConfig::frequency_bands[interp_high][0] && interp_high < MAX_GAIN_SUBBANDS - 1){
				++interp_high;
			}
			interp_low = interp_high - 1;
			float t = (freq - KinectConfig::frequency_bands[interp_low][0]) / (KinectConfig::frequency_bands[interp_high][0] - KinectConfig::frequency_bands[interp_low][0]);
			// special case 1.  Frequency is less than frequency_bands[0][0] - linear interpolate between 0 & frequency_bands[0][0]
			if (freq < KinectConfig::frequency_bands[0][0]){
				t = freq / KinectConfig::frequency_bands[0][0];
				for (int channel = 0; channel < MAX_MICROPHONES; ++channel){
					m_dynamic_gains[channel][index] = Utils::interpolate(zero, m_persistent_gains[channel][0], t);
				}
			}
			// special case 2, no need to interpolate
			else if (t == 0.f || t >= 1.f){
				int freq_index = t > 0.f ? interp_high : interp_low;
				for (int channel = 0; channel < MAX_MICROPHONES; ++channel){
					m_dynamic_gains[channel][index] = m_persistent_gains[channel][freq_index];
				}
			}
			// standard case  | here we need to interpolate the values
			else{
				for (int channel = 0; channel < MAX_MICROPHONES; ++channel){
					m_dynamic_gains[channel][index] = Utils::interpolate(m_persistent_gains[channel][interp_low], m_persistent_gains[channel][interp_high], t);
				}
			}
		}
	}
}
