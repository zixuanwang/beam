#include "Pipeline.h"

namespace Beam{
	Pipeline* Pipeline::p_instance = NULL;

	Pipeline::Pipeline() : m_noise_floor(20.0, 0.04, 30000.0, 0.0){
		// initialize noise suppressors.
		for (int channel = 0; channel < MAX_MICROPHONES; ++channel){
			m_ssl_noise_suppressor[channel].init(SAMPLE_RATE, FRAME_SIZE, 1.f, 10.f);
			m_pre_noise_suppressor[channel].init(SAMPLE_RATE, FRAME_SIZE, 1.f, 1.f);
			m_out_noise_suppressor[channel].init(SAMPLE_RATE, FRAME_SIZE, 1.f, 1.f);
		}
		DSPFilter::band_pass_mclt(m_ssl_band_pass_filter, 500.0 / SAMPLE_RATE, 1000.0 / SAMPLE_RATE, 2000.0 / SAMPLE_RATE, 3500.0 / SAMPLE_RATE);
		m_ssl.init(SAMPLE_RATE, FRAME_SIZE);
		// initialize m_frequency_filter.
		for (int sub = 0; sub < MAX_GAIN_SUBBANDS; ++sub){
			DSPFilter::band_pass_mclt(m_frequency_filter[sub], KinectConfig::frequency_bands[sub][0] / SAMPLE_RATE, KinectConfig::frequency_bands[sub][1] / SAMPLE_RATE, KinectConfig::frequency_bands[sub][1] / SAMPLE_RATE, KinectConfig::frequency_bands[sub][2] / SAMPLE_RATE);
		}
		// initialize m_working_frequency.
		m_working_frequency.assign(FRAME_SIZE, std::complex<float>(0.f, 0.f));
	}

	Pipeline* Pipeline::instance(){
		if (p_instance == NULL){
			p_instance = new Pipeline;
		}
		return p_instance;
	}

	void Pipeline::load_profile(){
		for (int beam = 0; beam < MAX_BEAMS; ++beam){
			for (int channel = 0; channel < MAX_MICROPHONES; ++channel){
				m_pcm_weights[beam][channel].assign(FRAME_SIZE, std::complex<float>(0.f, 0.f));
			}
		}
		// initialize kinect weights
		std::complex<float> zero(0.f, 0.f);
		float freq_step = (float)SAMPLE_RATE / FRAME_SIZE / 2.f;
		float freq_beg = freq_step / 2.f;
		for (int beam = 0; beam < MAX_BEAMS; ++beam){
			int interp_low = 0;
			int interp_high = 1;
			while (KinectConfig::kinect_descriptor.freq_low >= KinectConfig::kinect_weights.frequencies[interp_low] && interp_high < KinectConfig::kinect_weights.num_frequency_bins - 1){
				++interp_high;
				++interp_low;
			}
			for (int bin = 0; bin < FRAME_SIZE; ++bin){
				float freq = freq_beg + bin * freq_step;
				if (freq > KinectConfig::kinect_descriptor.freq_high){
					freq = KinectConfig::kinect_descriptor.freq_high;
				}
				while (freq >= KinectConfig::kinect_weights.frequencies[interp_high] && interp_high < KinectConfig::kinect_weights.num_frequency_bins - 1){
					++interp_high;
				}
				interp_low = interp_high - 1;
				float t = (freq - KinectConfig::kinect_weights.frequencies[interp_low]) / (KinectConfig::kinect_weights.frequencies[interp_high] - KinectConfig::kinect_weights.frequencies[interp_low]);
				// special case 1.  Frequency is less than dFreq_Lo --- set value to 0
				if (freq <= KinectConfig::kinect_descriptor.freq_low){
					for (int channel = 0; channel < MAX_MICROPHONES; ++channel){
						m_pcm_weights[beam][channel][bin] = zero;
					}
				}
				// Special Case 2 - interpolate between 0 and freqency
				else if (t < 0.f){
					int freq_index = KinectConfig::kinect_weights.frequencies[interp_low] > KinectConfig::kinect_descriptor.freq_low ? interp_low : interp_high;
					t = (freq - KinectConfig::kinect_descriptor.freq_low) / (KinectConfig::kinect_weights.frequencies[freq_index] - KinectConfig::kinect_descriptor.freq_low);
					for (int channel = 0; channel < MAX_MICROPHONES; ++channel){
						int weight_index = beam * FRAME_SIZE * MAX_MICROPHONES;
						weight_index += freq_index * KinectConfig::kinect_weights.num_channels;
						weight_index += channel;
						m_pcm_weights[beam][channel][bin] = Utils::interpolate(zero, KinectConfig::kinect_weights.weights[weight_index + KinectConfig::kinect_weights.num_channels], t);
					}
				}
				// special case 3, no need to interpolate
				else if (t == 0.f || t >= 1.f){
					int freq_index = t > 0.f ? interp_high : interp_low;
					for (int channel = 0; channel < MAX_MICROPHONES; ++channel){
						int weight_index = beam * KinectConfig::kinect_weights.num_frequency_bins * KinectConfig::kinect_weights.num_channels;
						weight_index += freq_index * KinectConfig::kinect_weights.num_channels;
						weight_index += channel;
						m_pcm_weights[beam][channel][bin] = KinectConfig::kinect_weights.weights[weight_index];
					}
				}
				// standard case  | here we need to interpolate the values
				else{
					for (int channel = 0; channel < MAX_MICROPHONES; ++channel){
						int weight_index = beam * KinectConfig::kinect_weights.num_frequency_bins * KinectConfig::kinect_weights.num_channels;
						weight_index += interp_low * KinectConfig::kinect_weights.num_channels;
						weight_index += channel;
						m_pcm_weights[beam][channel][bin] = Utils::interpolate(KinectConfig::kinect_weights.weights[weight_index], KinectConfig::kinect_weights.weights[weight_index + KinectConfig::kinect_weights.num_channels], t);
					}
				}
			}
		}
	}

	void Pipeline::preprocess(std::vector<std::complex<float> >* input){
		for (int channel = 0; channel < MAX_MICROPHONES; ++channel){
			m_pre_noise_suppressor[channel].phase_compensation(input[channel]);
		}
	}

	bool Pipeline::source_localize(std::vector<std::complex<float> >* input, double time, float* p_angle){
		//  Apply the SSL band pass filter to the input channels
		//  and have a separate copy of the input channels 
		//  for SSL purposes only
		for (int channel = 0; channel < MAX_MICROPHONES; ++channel){
			for (size_t bin = 0; bin < m_ssl_band_pass_filter.size(); ++bin){
				input[channel][bin] *= m_ssl_band_pass_filter[bin];
			}
		}
		//  Noise suppression
		//  We do heavy noise suppression as we don't care about the musical noises
		//  but we do cary to suppress stationaty noises
		double energy = 0.0;
		for (int channel = 0; channel < MAX_MICROPHONES; ++channel){
			m_ssl_noise_suppressor[channel].noise_compensation(input[channel]);
			energy += (double)Utils::computeRMS(input[channel]);
		}
		energy /= MAX_MICROPHONES;
		double floor = m_noise_floor.nextLevel(time, energy);
		if (energy > 5.290792 * floor){	// TODO: modify this threshold
			// sound signal
			float angle;
			float weight;
			float confidence;
			float std_dev;
			int valid;
			int num;
			m_ssl.process(input, &angle, &weight);
			std::cout << "weight: " << weight << std::endl;
			if (weight > 5e-6f){ // TODO: modify this threshold
				m_ssl.process_next_sample(time, angle, weight);
				m_ssl.get_average(time, p_angle, &confidence, &std_dev, &num, &valid);
				return true;
			}
		}
		return false;
	}

	void Pipeline::beamformer(std::vector<std::complex<float> >* input, std::vector<std::complex<float> >& output, double time){
		//  we have sound source detected - single beam mode
		if (m_confidence > SSL_BEAMCHANGE_CONFIDENCE_THRESHOLD){
			
		}
	}

	float Pipeline::smart_calibration(float sound_source, std::vector<std::complex<float> >* input){
		float channel_rms[MAX_MICROPHONES] = { 0.f };
		float est_channel_rms[MAX_MICROPHONES] = { 0.f };
		float est_gains[MAX_MICROPHONES] = { 0.f };
		float average_gain = 0.f;
		float sigma = -1.f;
		float theta_rad = 0.f;
		RCoords mic;
		//  Project microphones to the line pointing to the sound source
		//  Here we assume flat wave propagation from the sound source
		for (int channel = 0; channel < MAX_MICROPHONES; ++channel){
			Utils::c2r(mic, KinectConfig::kinect_descriptor.mic[channel].x, KinectConfig::kinect_descriptor.mic[channel].y, KinectConfig::kinect_descriptor.mic[channel].z);
			m_coordinates[channel] = mic.rho * cosf(sound_source - mic.fi) * cosf(mic.theta);
		}
		for (int sub = 0; sub < MAX_GAIN_SUBBANDS; ++sub){
			for (int channel = 0; channel < MAX_MICROPHONES; ++channel){
				for (int bin = 0; bin < FRAME_SIZE; ++bin){
					m_working_frequency[bin] = input[channel][bin] * m_frequency_filter[sub][bin];
				}
				channel_rms[channel] = Utils::computeRMS(m_working_frequency);
			}
			sigma = Utils::approx(m_coordinates, channel_rms, 1, m_coeff, MAX_MICROPHONES);
			float average_rms = 0.f;
			for (int channel = 0; channel < MAX_MICROPHONES; ++channel){
				average_rms += channel_rms[channel];
			}
			average_rms /= MAX_MICROPHONES;
			if (average_rms > FLT_MIN){
				sigma /= average_rms;
			}
			if (m_coeff[1] < 0.f){
				m_coeff[1] = 0.f;
				m_coeff[0] = average_rms;
				sigma = -1.f;
			}
			for (int channel = 0; channel < MAX_MICROPHONES; ++channel){
				est_channel_rms[channel] = m_coeff[1] * m_coordinates[channel] + m_coeff[0];
			}
			for (int channel = 0; channel < MAX_MICROPHONES; ++channel){
				est_gains[channel] = 1.f;
				if (channel_rms[channel] > FLT_MIN){
					est_gains[channel] = est_channel_rms[channel] / channel_rms[channel];
				}
			}
			for (int channel = 0; channel < MAX_MICROPHONES; ++channel){
				average_gain += est_gains[channel];
			}
			average_gain /= MAX_MICROPHONES;
			for (int channel = 0; channel < MAX_MICROPHONES; ++channel){
				est_gains[channel] /= average_gain;
			}
			float weight = 0.001f;
			for (int channel = 0; channel < MAX_MICROPHONES; ++channel){
				if (average_rms > FLT_MIN){
					weight = 0.001f * channel_rms[channel] / average_rms;
				}
				else{
					weight = 0.f;
				}
				float re = m_persistent_gains[channel][sub].real();
				float im = m_persistent_gains[channel][sub].imag();
				m_persistent_gains[channel][sub] = std::complex<float>(re + weight * (est_gains[channel] - re), im);
			}
		}
		return sigma;
	}

	void Pipeline::expand_gain(){
		std::vector<std::complex<float> > dynamic_gain[MAX_MICROPHONES];
		std::complex<float> zero(0.f, 0.f);

		// use MCLT
		float freq_step = (float)SAMPLE_RATE / FRAME_SIZE / 2.f;
		float freq_beg = freq_step / 2.f;
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
			// special case 1.  Frequency is less than g_FrequencyBands[0][0] - linear interpolate between 0 & g_FrequencyBands[0][0]
			if (freq < KinectConfig::frequency_bands[0][0]){
				t = freq / KinectConfig::frequency_bands[0][0];
				for (int channel = 0; channel < MAX_MICROPHONES; ++channel){
					std::complex<float> p1 = m_persistent_gains[channel][0];
					dynamic_gain[channel][index] = Utils::interpolate(zero, p1, t);
				}
			}
			// special case 2, no need to interpolate
			else if (t == 0.f || t >= 1.f){
				int freq_index = t > 0.f ? interp_high : interp_low;
				for (int channel = 0; channel < MAX_MICROPHONES; ++channel){
					dynamic_gain[channel][index] = m_persistent_gains[channel][freq_index];
				}
			}
			// standard case  | here we need to interpolate the values
			else{
				for (int channel = 0; channel < MAX_MICROPHONES; ++channel){
					dynamic_gain[channel][index] = Utils::interpolate(m_persistent_gains[channel][interp_low], m_persistent_gains[channel][interp_high], t);
				}
			}
		}
	}
}
