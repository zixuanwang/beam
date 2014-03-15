#include "SoundSourceLocalizer.h"

SoundSourceLocalizer::SoundSourceLocalizer(){

}

SoundSourceLocalizer::~SoundSourceLocalizer(){

}

void SoundSourceLocalizer::init(double sample_rate, int frame_size){
	m_sample_rate = sample_rate;
	m_frame_size = frame_size;
	std::complex<double> gain[MAX_MICROPHONES];
	double beg_angle = -1.0 * HALF_PI;
	double end_angle = HALF_PI;
	double step_angle = PI / (NUM_ANGLES - 1);
	for (int angle = 0; angle < NUM_ANGLES; ++angle){
		m_angle[angle] = beg_angle + step_angle * angle;
	}
	// synthetic data
	// computational bins
	m_start_bin = (int)floor(500.0 / m_sample_rate * 2.0 * m_frame_size + 0.5);
	m_end_bin = (int)floor(3500.0 / m_sample_rate * 2.0 * m_frame_size + 0.5);
	m_meas_bins = m_end_bin - m_start_bin + 1;
	for (int pair = 0; pair < (MAX_MICROPHONES - 1); ++pair){
		for (int angle = 0; angle < NUM_ANGLES; ++angle){
			m_delta[pair][angle] = std::unique_ptr<double[]>(new double[m_meas_bins]);
		}
	}
	for (int meas_bin = 0, bin = m_start_bin; meas_bin < m_meas_bins; ++meas_bin, ++bin){
		double freq = (double)bin * m_sample_rate / m_frame_size / 2.0;
		for (int angle = 0; angle < NUM_ANGLES; ++angle){
			double x = DISTANCE * cos(m_angle[angle]);
			double y = DISTANCE * sin(m_angle[angle]);
			double z = 0.0;
			for (int mic = 0; mic < MAX_MICROPHONES; ++mic){
				double dx = x - m_mic_array.mic_array[mic].x;
				double dy = y - m_mic_array.mic_array[mic].y;
				double dz = z - m_mic_array.mic_array[mic].z;
				double dist = sqrt(dx * dx + dy * dy + dz * dz);
				double xy_dist = sqrt(dx * dx + dy * dy);
				double gamma = atan2(dy, dx); // a bug in the source code on this line. L176
				double cappa = atan2(dz, xy_dist);
				gamma -= m_mic_array.mic_array[mic].direction * TO_RAD;
				cappa -= m_mic_array.mic_array[mic].elevation * TO_RAD;
				double cos_theta = cos(gamma) * cos(cappa);
				double im = -1.0 * TWO_PI * freq * dist / SOUND_SPEED;
				gain[mic] = std::complex<double>(cos(im), sin(im)) * Microphone::micRatio(cos_theta, freq) * (1.0 / dist);
			}
			for (int pair = 0; pair < (MAX_MICROPHONES - 1); ++pair){
				m_delta[pair][angle][meas_bin] = Math::normalize_angle(std::arg(gain[0]) - std::arg(gain[pair + 1]));
			}
		}
	}
}

void SoundSourceLocalizer::process(const std::vector<std::vector<std::complex<double> > >& input_channels, double* p_angle, double* p_weight){
	double delta[MAX_MICROPHONES - 1];
	double ssl_sum[NUM_ANGLES] = { 0.0 };
	int bin, meas_bin;
	for (bin = m_start_bin, meas_bin = 0; bin < m_end_bin; ++bin, ++meas_bin){
		double sample_amplitude = std::abs(input_channels[0][bin]);
		for (int pair = 0; pair < (MAX_MICROPHONES - 1); ++pair){
			delta[pair] = Math::normalize_angle(std::arg(input_channels[0][bin]) - std::arg(input_channels[pair + 1][bin]));
			sample_amplitude += std::abs(input_channels[pair + 1][bin]); // bug in the source code. L300
		}
		sample_amplitude /= (double)MAX_MICROPHONES;
		double min_dist = 1e10;
		int min_index = 0;
		for (int angle = 0; angle < NUM_ANGLES; ++angle){
			double dist = 0.0;
			for (int pair = 0; pair < (MAX_MICROPHONES - 1); ++pair){
				double c = Math::normalize_angle(delta[pair] - m_delta[pair][angle][bin]);
				dist += c * c;
			}
			if (dist < min_dist){
				min_dist = dist;
				min_index = angle;
			}
		}
		ssl_sum[min_index] += sample_amplitude;
	}
	auto iter = std::max_element(ssl_sum, ssl_sum + NUM_ANGLES);
	double weight_max = *iter;
	int max_index = (int)(iter - ssl_sum);
	// interpolate the distribution function maximum using second degree polynom.
	if ((max_index <= 0) || (max_index >= (NUM_ANGLES - 1))){
		*p_angle = m_angle[max_index];
	}
	else{
		std::vector<double> x = { m_angle[max_index - 1], m_angle[max_index], m_angle[max_index + 1] };
		std::vector<double> y = { ssl_sum[max_index - 1], ssl_sum[max_index], ssl_sum[max_index + 1] };
		*p_angle = Math::interpolateMax(x, y);
	}
	*p_weight = weight_max / m_meas_bins / NUM_ANGLES;
}