#ifndef SOUNDSOURCELOCALIZER_H_
#define SOUNDSOURCELOCALIZER_H_

#include <algorithm>
#include <memory>
#include <vector>
#include "Math.h"
#include "MicrophoneArray.h"


#define NUM_ANGLES 18
#define DISTANCE 1.5f

class SoundSourceLocalizer {
public:
	SoundSourceLocalizer();
	~SoundSourceLocalizer();
	void init(double sample_rate, int frame_size);
	void process(const std::vector<std::vector<std::complex<double> > >& input_channels, double* p_angle, double* p_weight);

private:
	// sythetic data.
	std::unique_ptr<double[]> m_delta[MAX_MICROPHONES - 1][NUM_ANGLES];
	double m_angle[NUM_ANGLES];
	int m_start_bin;
	int m_end_bin;
	int m_meas_bins;
	// common.
	MicrophoneArray m_mic_array;
	double m_sample_rate;
	int m_frame_size;
};

#endif /* SOUNDSOURCELOCALIZER_H_ */
