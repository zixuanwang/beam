#include "MVDRBeamformer.h"

namespace Beam{
	MVDRBeamformer::MVDRBeamformer(){
		m_nn.assign(FRAME_SIZE, arma::cx_fmat());
	}

	MVDRBeamformer::~MVDRBeamformer(){
	
	}

	void MVDRBeamformer::compute(std::vector<std::complex<float> >* input, std::vector<std::complex<float> >& output, float angle, float confidence, double time, bool voice){
		if (!voice){
			// noise frame, update noise covariance matrix
			for (int bin = 0; bin < FRAME_SIZE; ++bin){
				arma::cx_fmat nn = arma::zeros<arma::cx_fmat>(MAX_MICROPHONES, MAX_MICROPHONES);
				for (int i = 0; i < MAX_MICROPHONES; ++i){
					for (int j = 0; j < MAX_MICROPHONES; ++j){
						nn(i, j) = input[i][bin] * input[j][bin];
					}
				}
				if (m_nn[bin].n_elem == 0){
					m_nn[bin] = nn;
				}
				else{
					m_nn[bin] = m_nn[bin] * 0.9f + nn * 0.1f;
				}
			}
		}
		// compute time delay
		for (int bin = 0; bin < FRAME_SIZE; ++bin){
			float time_delay[MAX_MICROPHONES] = { 0.f };
			for (int channel = 0; channel < MAX_MICROPHONES; ++channel){
				float distance = KinectConfig::kinect_descriptor.mic[channel].y * sinf(angle);
				time_delay[channel] = distance / (float)SOUND_SPEED;
			}
			arma::cx_fmat d_h = arma::zeros<arma::cx_fmat>(MAX_MICROPHONES, 1);
			float rad_freq = (float)(-bin * TWO_PI * SAMPLE_RATE / FRAME_SIZE / 2.f);
			for (int i = 0; i < MAX_MICROPHONES; ++i){
				float angle = rad_freq * time_delay[i];
				d_h(i, 0) = arma::cx_float(cosf(angle), sinf(angle));
			}
			arma::cx_fmat d = d_h.t();
			std::cout << m_nn[bin] << std::endl;
			std::cout << arma::cond(m_nn[bin]) << std::endl;
			std::cout << d_h << std::endl;
			arma::cx_fmat nn_inv_d_h = arma::solve(m_nn[bin], d_h);
			arma::cx_fmat denom_mat = d * nn_inv_d_h;
			float denom = denom_mat(0, 0).real();
			std::complex<float> sum(0.f, 0.f);
			for (int channel = 0; channel < MAX_MICROPHONES; ++channel){
				sum += input[channel][bin] * nn_inv_d_h(channel, 0);
			}
			output[bin] = sum / denom;
		}
	}
}