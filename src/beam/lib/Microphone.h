#ifndef MICROPHONE_H_
#define MICROPHONE_H_

#include <complex>
#include "Math.h"

#define SOUND_SPEED 342.0f;

class Microphone {
public:
	Microphone();
	~Microphone();
	/// returns the complex gain for omni microphones.
	static std::complex<float> micRatio(float cos_theta, float freq);
	/// returns the complex gain for ideal microphone given alpha and beta
	static std::complex<float> microphoneDirectivity(float freq, float cos_theta, float alpha, float beta);
	/// returns the complex gain for ideal gradient microphone.
	static std::complex<float> gradientMicrophoneDirectivity(float freq, float cos_theta);


	/// coordinates of the microphone.
	float x;
	float y;
	float z;
	/// direction of the microphone.
	float direction;
	/// elevation of the microphone.
	float elevation;
};

#endif /* MICROPHONE_H_ */
