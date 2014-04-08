#ifndef SOUNDSOURCELOCALIZER_H_
#define SOUNDSOURCELOCALIZER_H_

#include <algorithm>
#include <list>
#include <memory>
#include "KinectConfig.h"

namespace Beam{
#define NUM_ANGLES 18
#define DISTANCE 1.5f
#define MAX_COORD_SAMPLES 40
#define NUM_CLUSTERS 36
#define SSL_CONFIDENT_MEASUREMENTS 1.5f
#define SSL_MEASUREMENT_DEVIATION 0.032656f
#define	SSL_MEASUREMENT_LIFETIME 2.5
#define SSL_CONTRAST_THRESHOLD 0.2716f
#define	SSL_BEAMCHANGE_CONFIDENCE_THRESHOLD 0.431948f
	class SoundSourceLocalizer {
	public:
		SoundSourceLocalizer();
		~SoundSourceLocalizer();
		void init(float sample_rate, int frame_size);
		void process(std::vector<std::complex<float> >* input, float* p_angle, float* p_weight);
		void process_next_sample(double time, float next_point, float weight);
		/// filtering the angle.
		void get_average(double time, float* p_average, float* p_confidence, float* p_std_dev, int* p_num, int* p_valid);
	private:
		struct CoordsSample{
			double time;
			float point;
			float weight;
		};
		struct Cluster{
			float average;
			float std_dev;
			float weight;
			float points[MAX_COORD_SAMPLES];
			float weights[MAX_COORD_SAMPLES];
			float diffs[MAX_COORD_SAMPLES];
			int num_points;
			int valid_points;
		};
		// sythetic data.
		std::unique_ptr<float[]> m_delta[MAX_MICROPHONES - 1][NUM_ANGLES];
		float m_angle[NUM_ANGLES];
		int m_start_bin;
		int m_end_bin;
		int m_meas_bins;
		// common.
		float m_sample_rate;
		int m_frame_size;
		// record samples
		std::list<Beam::SoundSourceLocalizer::CoordsSample> m_coord_samples;
		bool m_new_sample;
		double m_last_time; // last time the angle is known.
		Beam::SoundSourceLocalizer::Cluster m_sample_cluster[NUM_CLUSTERS];
		float m_upper_boundary[NUM_CLUSTERS];
		float m_lower_boundary[NUM_CLUSTERS];
	};
}

#endif /* SOUNDSOURCELOCALIZER_H_ */
