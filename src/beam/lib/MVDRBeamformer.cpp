#include "MVDRBeamformer.h"

namespace Beam {
MVDRBeamformer::MVDRBeamformer() {
	m_nn.assign(FRAME_SIZE, Eigen::Matrix<std::complex<float>, MAX_MICROPHONES, MAX_MICROPHONES>::Identity());
}

MVDRBeamformer::~MVDRBeamformer() {

}

void MVDRBeamformer::compute(std::vector<std::complex<float> >* input,
		std::vector<std::complex<float> >& output, float angle,
		float confidence, double time, bool voice) {
	if (!voice) {
		// noise frame, update noise covariance matrix
		for (int bin = 0; bin < FRAME_SIZE; ++bin) {
			Eigen::Matrix<std::complex<float>, MAX_MICROPHONES, 1> n;
			for (int i = 0; i < MAX_MICROPHONES; ++i) {
				n(i, 0) = input[i][bin];
			}
			Eigen::Matrix<std::complex<float>, MAX_MICROPHONES, MAX_MICROPHONES> nn = n * n.adjoint();
			m_nn[bin] = m_nn[bin] * 0.99f + nn * 0.01f;
		}
	}
	// compute time delay
	for (int bin = 0; bin < FRAME_SIZE; ++bin) {
		Eigen::Matrix<std::complex<float>, MAX_MICROPHONES, MAX_MICROPHONES> inverse;
		std::complex<float> determinant;
		bool invertible;
		m_nn[bin].computeInverseAndDetWithCheck(inverse,determinant,invertible);
		if (invertible) {
			float time_delay[MAX_MICROPHONES] = { 0.f };
			for (int channel = 0; channel < MAX_MICROPHONES; ++channel) {
				float distance = KinectConfig::kinect_descriptor.mic[channel].y
						* sinf(angle);
				time_delay[channel] = distance / (float) SOUND_SPEED;
			}
			Eigen::Matrix<std::complex<float>, MAX_MICROPHONES, 1> d_h;
			float rad_freq = (float) (-bin * TWO_PI * SAMPLE_RATE / FRAME_SIZE
					/ 2.f);
			for (int i = 0; i < MAX_MICROPHONES; ++i) {
				float angle = rad_freq * time_delay[i];
				d_h(i, 0) = std::complex<float>(cosf(angle), sinf(angle));
			}
			Eigen::Matrix<std::complex<float>, 1, MAX_MICROPHONES> d = d_h.adjoint();
			Eigen::Matrix<std::complex<float>, MAX_MICROPHONES, 1> nn_inv_d_h = inverse * d_h;
			Eigen::Matrix<std::complex<float>, 1, 1> denom_mat = d * nn_inv_d_h;
			float denom = denom_mat(0, 0).real();
			std::complex<float> sum(0.f, 0.f);
			for (int channel = 0; channel < MAX_MICROPHONES; ++channel) {
				sum += input[channel][bin] * nn_inv_d_h(channel, 0);
			}
			output[bin] = sum / denom;
		}
	}
}
}
