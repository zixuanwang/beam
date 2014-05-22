#include "SoundSourceLocalizer.h"

#include <iostream>

namespace Beam{
	SoundSourceLocalizer::SoundSourceLocalizer(){
		m_new_sample = false;
		m_last_time = 0.0;
		float step = (float)TWO_PI / NUM_CLUSTERS;
		m_lower_boundary[0] = (float)-PI;
		m_upper_boundary[0] = m_lower_boundary[0] + 2.f * step;
		for (int i = 1; i < NUM_CLUSTERS; ++i){
			m_lower_boundary[i] = m_lower_boundary[i - 1] + step;
			m_upper_boundary[i] = m_upper_boundary[i - 1] + step;
		}
		m_confidence = 0.f;
		m_std_dev = 0.f;
		m_average = 0.f;
		m_valid = 0;
		m_num = 0;
	}

	SoundSourceLocalizer::~SoundSourceLocalizer(){

	}

	void SoundSourceLocalizer::init(float sample_rate, int frame_size){
		m_sample_rate = sample_rate;
		m_frame_size = frame_size;
		std::complex<float> gain[MAX_MICROPHONES];
		double beg_angle = -HALF_PI;
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
					float dx = x - KinectConfig::kinect_descriptor.mic[mic].x;
					float dy = y - KinectConfig::kinect_descriptor.mic[mic].y;
					float dz = z - KinectConfig::kinect_descriptor.mic[mic].z;
					float dist = sqrtf(dx * dx + dy * dy + dz * dz);
					float xy_dist = sqrtf(dx * dx + dy * dy);
					float gamma = atan2(dy, dx);
					float cappa = atan2(dz, xy_dist);
					gamma -= (float)(KinectConfig::kinect_descriptor.mic[mic].direction * TO_RAD);
					cappa -= (float)(KinectConfig::kinect_descriptor.mic[mic].elevation * TO_RAD);
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
			float sample_amplitude = Utils::abs_complex(input[0][bin]);
			for (int pair = 0; pair < (MAX_MICROPHONES - 1); ++pair){
				delta[pair] = Utils::normalize_angle(std::arg(input[0][bin]) - std::arg(input[pair + 1][bin]));
				sample_amplitude += Utils::abs_complex(input[pair + 1][bin]); // bug in the source code. L300
			}
			sample_amplitude /= (float)MAX_MICROPHONES;
			float min_dist = 1e10f;
			int min_index = 0;
			for (int angle = 0; angle < NUM_ANGLES; ++angle){
				float dist = 0.0;
				for (int pair = 0; pair < (MAX_MICROPHONES - 1); ++pair){
					float c = Utils::normalize_angle(delta[pair] - m_delta[pair][angle][meas_bin]);
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

	void SoundSourceLocalizer::process_next_sample(double time, float next_point, float weight){
		// add next point to the measurements queue
		// TODO check the limits.
		//Utils::limit(weight, 0.f, 5.f);
		m_coord_samples.push_back(CoordsSample{ time, next_point, weight });
		if (m_coord_samples.size() > MAX_COORD_SAMPLES){
			m_coord_samples.pop_front();
		}
		m_new_sample = true;
	}

	void SoundSourceLocalizer::get_average(double time, float* p_average, float* p_confidence, float* p_std_dev, int* p_num, int* p_valid){
		*p_average = 0.f;
		*p_confidence = 0.f;
		*p_std_dev = 0.f;
		*p_num = 0;
		*p_valid = 0;
		if (!m_new_sample && (time - m_last_time < 0.2)){
			// no new measurements and fresh data.
			*p_average = m_average;
			*p_confidence = m_confidence;
			*p_std_dev = m_std_dev;
			*p_num = m_num;
			*p_valid = m_valid;
			return;
		}
		m_last_time = time;
		m_new_sample = false;
		//  cleanup the measurements array - remove meassurements older than m_life_time
		double last_valid_time = time - SSL_MEASUREMENT_LIFETIME;
		//  find the begin of the valid measurement.
		auto valid_iter = m_coord_samples.begin();
		for (; valid_iter != m_coord_samples.end(); ++valid_iter){
			if (valid_iter->time >= last_valid_time) break;
		}
		if (m_coord_samples.empty() || valid_iter == m_coord_samples.end()){
			// queue is empty
			return;
		}
		double last_measurement_time = m_coord_samples.back().time;
		//  prepare the clustering
		for (int i = 0; i < NUM_CLUSTERS; ++i){
			m_sample_cluster[i].average = 0.f;
			m_sample_cluster[i].weight = 0.f;
			m_sample_cluster[i].std_dev = 0.f;
			m_sample_cluster[i].num_points = 0;
		}
		//  cluster the measurements - horizontal plane only
		for (auto iter = valid_iter; iter != m_coord_samples.end(); ++iter){
			float angle = iter->point;
			for (int i = 0; i < NUM_CLUSTERS; ++i){
				if (angle < m_upper_boundary[i] && angle >= m_lower_boundary[i]){
					int index = m_sample_cluster[i].num_points;
					m_sample_cluster[i].points[index] = angle;
					m_sample_cluster[i].weights[index] = iter->weight;
					++m_sample_cluster[i].num_points;
				}
			}
		}
		//  process the measurements cluster by cluster
		for (int i = 0; i < NUM_CLUSTERS; ++i){
			for (int j = 0; j < m_sample_cluster[i].num_points; ++j){
				m_sample_cluster[i].average += m_sample_cluster[i].points[j] * m_sample_cluster[i].weights[j];
				m_sample_cluster[i].weight += m_sample_cluster[i].weights[j];
			}
			if (m_sample_cluster[i].num_points > 0){
				m_sample_cluster[i].average /= m_sample_cluster[i].weight;
			}
		}
		//  Calculate the standard deviation, cluster by cluster
		for (int i = 0; i < NUM_CLUSTERS; ++i){
			for (int j = 0; j < m_sample_cluster[i].num_points; ++j){
				float diff = fabs(m_sample_cluster[i].points[j] - m_sample_cluster[i].average);
				m_sample_cluster[i].diffs[j] = diff;
				m_sample_cluster[i].std_dev += diff * diff * m_sample_cluster[i].weights[j];
			}
			if (m_sample_cluster[i].num_points > 0){
				m_sample_cluster[i].std_dev = sqrtf(m_sample_cluster[i].std_dev / m_sample_cluster[i].weight);
			}
		}
		//  recalculate the average without measurements out of +/- 2.0 StdDev
		//  makes sense only if we have enough measurements
		for (int i = 0; i < NUM_CLUSTERS; ++i){
			if (m_sample_cluster[i].num_points > 10){
				float window = 2.f * m_sample_cluster[i].std_dev;
				float temp_avg = 0.f;
				float temp_weight = 0.f;
				m_sample_cluster[i].valid_points = 0;
				//  Recalculate the average
				for (int j = 0; j < m_sample_cluster[i].num_points; ++j){
					if (m_sample_cluster[i].diffs[j] <= window){
						temp_avg += m_sample_cluster[i].points[j] * m_sample_cluster[i].weights[j];
						temp_weight += m_sample_cluster[i].weights[j];
						++m_sample_cluster[i].valid_points;
					}
				}
				//  Recalculate the standard deviation
				float temp_std_dev = 0.f;
				if (m_sample_cluster[i].valid_points > 0){
					temp_avg /= temp_weight;
					for (int j = 0; j < m_sample_cluster[i].num_points; ++j){
						if (m_sample_cluster[i].diffs[j] <= window){
							float diff = fabs(m_sample_cluster[i].points[j] - temp_avg);
							m_sample_cluster[i].diffs[j] = diff;
							temp_std_dev += diff * diff * m_sample_cluster[i].weights[j];
						}
					}
				}
				//  Use the new average only if we didn't remove more than
				//  one half of the measurements - in the other case 
				//  something is very wrong with our measurements!
				//  Bad (or very different than Gausian) distribution
				//  model (two sound sources, for example). Theorethically 
				//  we should use > 95% of the measurements
				if ((m_sample_cluster[i].valid_points > m_sample_cluster[i].num_points / 2) &&
					(m_sample_cluster[i].valid_points > 1)){
					temp_std_dev = sqrtf(temp_std_dev / temp_weight);
					m_sample_cluster[i].average = temp_avg;
					m_sample_cluster[i].std_dev = temp_std_dev;
				}
				else{
					m_sample_cluster[i].valid_points = m_sample_cluster[i].num_points;
				}
			}
			else{
				m_sample_cluster[i].valid_points = m_sample_cluster[i].num_points;
			}
		}
		float weight = m_sample_cluster[0].weight;
		float average = m_sample_cluster[0].average;
		m_std_dev = m_sample_cluster[0].std_dev;
		m_valid = m_sample_cluster[0].valid_points;
		m_num = m_sample_cluster[0].num_points;
		for (int i = 1; i < NUM_CLUSTERS; ++i){
			if (m_sample_cluster[i].weight > weight){
				weight = m_sample_cluster[i].weight;
				average = m_sample_cluster[i].average;
				m_std_dev = m_sample_cluster[i].std_dev;
				m_valid = m_sample_cluster[i].valid_points;
				m_num = m_sample_cluster[i].num_points;
			}
		}
		m_average = average;
		//  claculate the confidence level:
		//  it depends on:
		//      the number of measurements [10 - 100%], 
		//      the standard deviation [0.2 rad - 100%] and 
		//      the time of the last measurement [ m_dLifeTime / 2 - 100%]
		float confidence1, confidence2, confidence3;
		if ((float)m_valid >= SSL_CONFIDENT_MEASUREMENTS){
			confidence1 = 1.f;
		}
		else{
			confidence1 = (float)m_valid / SSL_CONFIDENT_MEASUREMENTS;
		}
		if (m_std_dev <= SSL_MEASUREMENT_DEVIATION){
			confidence2 = 1.f;
		}
		else{
			confidence2 = SSL_MEASUREMENT_DEVIATION / m_std_dev;
		}
		if (confidence2 > 1.f) confidence2 = 1.f;
		confidence3 = (float)((SSL_MEASUREMENT_LIFETIME - (time - last_measurement_time)) / SSL_MEASUREMENT_LIFETIME * 2.0);
		if (confidence3 > 1.f){
			confidence3 = 1.f;
		}
		m_confidence = confidence1 * confidence2 * confidence3;
		Utils::limit(m_confidence, 0.f, 1.f);
		*p_average = m_average;
		*p_confidence = m_confidence;
		*p_std_dev = m_std_dev;
		*p_num = m_num;
		*p_valid = m_valid;
	}
}