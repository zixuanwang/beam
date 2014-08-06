#include "GSCBeamformer.h"

namespace Beam{
	GSCBeamformer::GSCBeamformer(){
		for (int i = 0; i < MAX_MICROPHONES; ++i){
			m_bm[i].assign(FRAME_SIZE, std::complex<float>(1.f, 0.f));
			m_y[i].assign(FRAME_SIZE, std::complex<float>(0.f, 0.f));
			m_w[i].assign(FRAME_SIZE, std::complex<float>(1.f / MAX_MICROPHONES, 0.f));
		}
	}

	GSCBeamformer::~GSCBeamformer(){
	
	}

	void GSCBeamformer::compute(std::vector<std::complex<float> >* input, std::vector<std::complex<float> >& output, float angle, float confidence, double time, bool voice){
		// compute time delay
		std::vector<std::complex<float> > fixed_beamformer(FRAME_SIZE);
		float time_delay[MAX_MICROPHONES] = { 0.f };
		for (int channel = 0; channel < MAX_MICROPHONES; ++channel){
			float distance = KinectConfig::kinect_descriptor.mic[channel].y * sinf(angle);
			time_delay[channel] = distance / (float)SOUND_SPEED;
		}
		for (int bin = 0; bin < FRAME_SIZE; ++bin){
			std::complex<float> sum(0.f, 0.f);
			float rad_freq = (float)(-bin * TWO_PI * SAMPLE_RATE / FRAME_SIZE / 2.f);
			for (int channel = 0; channel < MAX_MICROPHONES; ++channel){
				float v = (float)(rad_freq * time_delay[channel]);
				sum += input[channel][bin] * std::complex<float>(cosf(v), sinf(v));
			}
			sum /= (float)MAX_MICROPHONES;
			fixed_beamformer[bin] = sum;
		}
		if (voice){
			// update adaptive blocking matrix
			for (int i = 0; i < MAX_MICROPHONES; ++i){
				for (int bin = 0; bin < FRAME_SIZE; ++bin){
					float norm = std::norm(fixed_beamformer[bin]);
					if (norm != 0.f){
						m_bm[i][bin] += (input[i][bin] - m_bm[i][bin] * fixed_beamformer[bin]) * fixed_beamformer[bin] / norm;
					}
				}
			}
		}
		for (int i = 0; i < MAX_MICROPHONES; ++i){
			for (int bin = 0; bin < FRAME_SIZE; ++bin){
				m_y[i][bin] = input[i][bin] - m_bm[i][bin] * fixed_beamformer[bin];
			}
		}
		for (int bin = 0; bin < FRAME_SIZE; ++bin){
			output[bin] = m_y[0][bin];
		}
		/*
		if (!voice){
			// update adaptive input canceller
			for (int bin = 0; bin < FRAME_SIZE; ++bin){
				std::complex<float> sum(0.f, 0.f);
				std::complex<float> norm(0.f, 0.f);
				for (int i = 0; i < MAX_MICROPHONES; ++i){
					sum += m_w[i][bin] * m_y[i][bin];
					norm += std::norm(m_y[i][bin]);
				}
				for (int i = 0; i < MAX_MICROPHONES; ++i){
					m_w[i][bin] -= (fixed_beamformer[bin] - sum) * m_y[i][bin] / norm;
				}
			}
		}
		for (int bin = 0; bin < FRAME_SIZE; ++bin){
			std::complex<float> sum(0.f, 0.f);
			for (int i = 0; i < MAX_MICROPHONES; ++i){
				sum += m_w[i][bin] * m_y[i][bin];
			}
			output[bin] = fixed_beamformer[bin] - sum;
		}
		*/
	}
}