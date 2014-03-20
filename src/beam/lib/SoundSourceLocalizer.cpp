#include "SoundSourceLocalizer.h"

namespace Beam{
	SoundSourceLocalizer::SoundSourceLocalizer(){

	}

	SoundSourceLocalizer::~SoundSourceLocalizer(){

	}

	void SoundSourceLocalizer::init(float sample_rate, int frame_size){
		m_sample_rate = sample_rate;
		m_frame_size = frame_size;
		std::complex<float> gain[MAX_MICROPHONES];
		double beg_angle = -1.0 * HALF_PI;
		double end_angle = HALF_PI;
		double step_angle = PI / (NUM_ANGLES - 1);
		for (int angle = 0; angle < NUM_ANGLES; ++angle){
			m_angle[angle] = (float)(beg_angle + step_angle * angle);
		}
		// synthetic data
		// computational bins
		m_start_bin = (int)floorf(500.f / m_sample_rate * 2.f * m_frame_size + 0.5f);
		m_end_bin = (int)floorf(3500.f / m_sample_rate * 2.f * m_frame_size + 0.5f);
		m_meas_bins = m_end_bin - m_start_bin + 1;
		for (int pair = 0; pair < (MAX_MICROPHONES - 1); ++pair){
			for (int angle = 0; angle < NUM_ANGLES; ++angle){
				m_delta[pair][angle] = std::unique_ptr<float[]>(new float[m_meas_bins]);
			}
		}
		for (int meas_bin = 0, bin = m_start_bin; meas_bin < m_meas_bins; ++meas_bin, ++bin){
			float freq = (float)bin * m_sample_rate / m_frame_size / 2.f;
			for (int angle = 0; angle < NUM_ANGLES; ++angle){
				float x = DISTANCE * cosf(m_angle[angle]);
				float y = DISTANCE * sinf(m_angle[angle]);
				float z = 0.f;
				for (int mic = 0; mic < MAX_MICROPHONES; ++mic){
					float dx = x - m_mic_array.mic_array[mic].x;
					float dy = y - m_mic_array.mic_array[mic].y;
					float dz = z - m_mic_array.mic_array[mic].z;
					float dist = sqrtf(dx * dx + dy * dy + dz * dz);
					float xy_dist = sqrtf(dx * dx + dy * dy);
					float gamma = atan2(dy, dx);
					float cappa = atan2(dz, xy_dist);
					gamma -= (float)(m_mic_array.mic_array[mic].direction * TO_RAD);
					cappa -= (float)(m_mic_array.mic_array[mic].elevation * TO_RAD);
					float cos_theta = cosf(gamma) * cosf(cappa);
					float im = (float)(-TWO_PI * freq * dist / SOUND_SPEED);
					gain[mic] = std::complex<float>(cosf(im), sinf(im)) * Microphone::micRatio(cos_theta, freq) * (1.f / dist);
				}
				for (int pair = 0; pair < (MAX_MICROPHONES - 1); ++pair){
					m_delta[pair][angle][meas_bin] = Utils::normalize_angle(std::arg(gain[0]) - std::arg(gain[pair + 1]));
				}
			}
		}
	}

	void SoundSourceLocalizer::process(std::vector<std::complex<float> >* input, float* p_angle, float* p_weight){
		float delta[MAX_MICROPHONES - 1];
		float ssl_sum[NUM_ANGLES] = { 0.f };
		int bin, meas_bin;
		for (bin = m_start_bin, meas_bin = 0; bin < m_end_bin; ++bin, ++meas_bin){
			float sample_amplitude = std::abs(input[0][bin]);
			for (int pair = 0; pair < (MAX_MICROPHONES - 1); ++pair){
				delta[pair] = Utils::normalize_angle(std::arg(input[0][bin]) - std::arg(input[pair + 1][bin]));
				sample_amplitude += std::abs(input[pair + 1][bin]); // bug in the source code. L300
			}
			sample_amplitude /= (float)MAX_MICROPHONES;
			float min_dist = 1e10f;
			int min_index = 0;
			for (int angle = 0; angle < NUM_ANGLES; ++angle){
				float dist = 0.0;
				for (int pair = 0; pair < (MAX_MICROPHONES - 1); ++pair){
					float c = Utils::normalize_angle(delta[pair] - m_delta[pair][angle][bin]);
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
		float weight_max = *iter;
		int max_index = (int)(iter - ssl_sum);
		// interpolate the distribution function maximum using second degree polynom.
		if ((max_index <= 0) || (max_index >= (NUM_ANGLES - 1))){
			*p_angle = m_angle[max_index];
		}
		else{
			std::vector<float> x = { m_angle[max_index - 1], m_angle[max_index], m_angle[max_index + 1] };
			std::vector<float> y = { ssl_sum[max_index - 1], ssl_sum[max_index], ssl_sum[max_index + 1] };
			*p_angle = Utils::interpolate_max(x, y);
		}
		*p_weight = weight_max / m_meas_bins / NUM_ANGLES;
	}
}
