#include "Pipeline.h"

namespace Beam{
	Pipeline* Pipeline::p_instance = NULL;

	Pipeline::Pipeline() : m_noise_floor(20.0, 0.04, 30000.0, 0.0){
		for (int channel = 0; channel < MAX_MICROPHONES; ++channel){
			m_ssl_noise_suppressor[channel].init(SAMPLE_RATE, FRAME_SIZE, 1.f, 10.f);
			m_pre_noise_suppressor[channel].init(SAMPLE_RATE, FRAME_SIZE, 1.f, 1.f);
			m_out_noise_suppressor[channel].init(SAMPLE_RATE, FRAME_SIZE, 1.f, 1.f);
		}
		DSPFilter::band_pass_mclt(m_ssl_band_pass_filter, 500.0 / SAMPLE_RATE, 1000.0 / SAMPLE_RATE, 2000.0 / SAMPLE_RATE, 3500.0 / SAMPLE_RATE);
		m_ssl.init(SAMPLE_RATE, FRAME_SIZE);
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
		for (int beam = 0; beam < KinectConfig::kinect_weights.num_beams; ++beam){
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
				if (freq <= KinectConfig::kinect_descriptor.freq_low){
					for (int channel = 0; channel < MAX_MICROPHONES; ++channel){
						m_pcm_weights[beam][channel][bin] *= 0.f;
					}
				}
				else if (t < 0.f){
					int freq_index = KinectConfig::kinect_weights.frequencies[interp_low] > KinectConfig::kinect_descriptor.freq_low ? interp_low : interp_high;
					t = (freq - KinectConfig::kinect_descriptor.freq_low) / (KinectConfig::kinect_weights.frequencies[freq_index] - KinectConfig::kinect_descriptor.freq_low);
					for (int channel = 0; channel < MAX_MICROPHONES; ++channel){
						int weight_index = beam * KinectConfig::kinect_weights.num_frequency_bins * KinectConfig::kinect_weights.num_channels;
						weight_index += freq_index * KinectConfig::kinect_weights.num_channels;
						weight_index += channel;
						m_pcm_weights[beam][channel][bin] = Utils::interpolate(zero, KinectConfig::kinect_weights.weights[weight_index + KinectConfig::kinect_weights.num_channels], t);
					}
				}
				else if (t == 0.f || t >= 1.f){
					int freq_index = t > 0.f ? interp_high : interp_low;
					for (int channel = 0; channel < MAX_MICROPHONES; ++channel){
						int weight_index = beam * KinectConfig::kinect_weights.num_frequency_bins * KinectConfig::kinect_weights.num_channels;
						weight_index += freq_index * KinectConfig::kinect_weights.num_channels;
						weight_index += channel;
						m_pcm_weights[beam][channel][bin] = KinectConfig::kinect_weights.weights[weight_index];
					}
				}
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
}
