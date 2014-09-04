#include "DelaySumBeamformer.h"

namespace Beam{
	DelaySumBeamformer::DelaySumBeamformer(){

	}

	DelaySumBeamformer::~DelaySumBeamformer(){
	
	}

	void DelaySumBeamformer::compute(std::vector<std::complex<float> >* input, std::vector<std::complex<float> >& output, float angle, float confidence, double time){
		// compute time delay
		float time_delay[AVALABLE_MICROPHONES] = { 0.f };
		for (int channel = 0; channel < AVALABLE_MICROPHONES; ++channel){
			float distance = KinectConfig::kinect_descriptor.mic[channel].y * sinf(angle);
			time_delay[channel] = distance / (float)SOUND_SPEED;
		}
		for (int bin = 0; bin < FRAME_SIZE; ++bin){
			std::complex<float> sum(0.f, 0.f);
			float rad_freq = (float)(-bin * TWO_PI * SAMPLE_RATE / FRAME_SIZE / 2.f);
			for (int channel = 0; channel < AVALABLE_MICROPHONES; ++channel){
				float v = (float)(rad_freq * time_delay[channel]);
				sum += input[channel][bin] * std::complex<float>(cosf(v), sinf(v));
			}
			sum /= (float)AVALABLE_MICROPHONES;
			output[bin] = sum;
		}
	}
}