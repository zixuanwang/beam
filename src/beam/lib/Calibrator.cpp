#include "Calibrator.h"

namespace Beam{
	Calibrator::Calibrator(){
		// initialize m_frequency_filter.
		for (int sub = 0; sub < MAX_GAIN_SUBBANDS; ++sub){
			DSPFilter::band_pass_mclt(m_frequency_filter[sub], KinectConfig::frequency_bands[sub][1] / SAMPLE_RATE, KinectConfig::frequency_bands[sub][0] / SAMPLE_RATE, KinectConfig::frequency_bands[sub][0] / SAMPLE_RATE, KinectConfig::frequency_bands[sub][2] / SAMPLE_RATE);
		}
		// initialize m_working_frequency.
		m_working_frequency.assign(FRAME_SIZE, std::complex<float>(0.f, 0.f));
	}

	Calibrator::~Calibrator(){
	
	}

	float Calibrator::calibrate(float sound_source, std::vector<std::complex<float> >* input, std::complex<float> persistent_gains[MAX_MICROPHONES][MAX_GAIN_SUBBANDS]){
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
			// TODO: maybe a bug here 
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
			average_gain = 0.f;
			for (int channel = 0; channel < MAX_MICROPHONES; ++channel){
				average_gain += est_gains[channel];
			}
			average_gain /= MAX_MICROPHONES;
			for (int channel = 0; channel < MAX_MICROPHONES; ++channel){
				est_gains[channel] /= average_gain;
			}
			for (int channel = 0; channel < MAX_MICROPHONES; ++channel){
				if (!((est_gains[channel] > 0.5f) && (est_gains[channel] < 2.f)))
					return -1.f;
			}
			float weight = 0.001f;
			for (int channel = 0; channel < MAX_MICROPHONES; ++channel){
				if (average_rms > FLT_MIN){
					weight = 0.001f * channel_rms[channel] / average_rms;
				}
				else{
					weight = 0.f;
				}
				float re = persistent_gains[channel][sub].real();
				float im = persistent_gains[channel][sub].imag();
				persistent_gains[channel][sub].real(re + weight * (est_gains[channel] - re));
			}
		}
		return sigma;
	}
}