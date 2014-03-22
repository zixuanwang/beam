#include "Pipeline.h"

namespace Beam{
	Pipeline* Pipeline::p_instance = NULL;

	Pipeline::Pipeline() : m_noise_floor(20.0, 0.04, 30000.0, 0.0){
		for (int channel = 0; channel < MAX_MICROPHONES; ++channel){
			m_ssl_noise_suppressor[channel].init(SAMPLE_RATE, FRAME_SIZE, 1.f, 10.f);
			m_pre_noise_suppressor[channel].init(SAMPLE_RATE, FRAME_SIZE, 1.f, 1.f);
		}
		DSPFilter::band_pass(m_ssl_band_pass_filter, 500.0 / SAMPLE_RATE, 1000.0 / SAMPLE_RATE, 2000.0 / SAMPLE_RATE, 3500.0 / SAMPLE_RATE);
	}

	Pipeline* Pipeline::instance(){
		if (p_instance == NULL){
			p_instance = new Pipeline;
		}
		return p_instance;
	}

	void Pipeline::load_profile(){

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
		// TODO: change the hard coding parameters.
		if ((energy > 5.290792 * floor) && (energy > 0.005)){
			// sound signal
			SoundSourceLocalizer ssl;
			ssl.init(SAMPLE_RATE, FRAME_SIZE);
			float weight;
			ssl.process(input, p_angle, &weight);
			return true;
		}
		return false;
	}
}
