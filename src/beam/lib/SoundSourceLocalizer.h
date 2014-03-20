#ifndef SOUNDSOURCELOCALIZER_H_
#define SOUNDSOURCELOCALIZER_H_

#include <algorithm>
#include <memory>
#include "MicrophoneArray.h"

namespace Beam{
#define NUM_ANGLES 18
#define DISTANCE 1.5f
	class SoundSourceLocalizer {
	public:
		SoundSourceLocalizer();
		~SoundSourceLocalizer();
		void init(float sample_rate, int frame_size);
		void process(std::vector<std::complex<float> >* input, float* p_angle, float* p_weight);

	private:
		// sythetic data.
		std::unique_ptr<float[]> m_delta[MAX_MICROPHONES - 1][NUM_ANGLES];
		float m_angle[NUM_ANGLES];
		int m_start_bin;
		int m_end_bin;
		int m_meas_bins;
		// common.
		MicrophoneArray m_mic_array;
		float m_sample_rate;
		int m_frame_size;
	};
}

#endif /* SOUNDSOURCELOCALIZER_H_ */
