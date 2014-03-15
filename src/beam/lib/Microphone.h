#ifndef MICROPHONE_H_
#define MICROPHONE_H_

#include <complex>
#include "Math.h"

#define SOUND_SPEED 342.0;

class Microphone {
public:
	Microphone();
	~Microphone();
	/// returns the complex gain for omni microphones.
	static std::complex<double> micRatio(double cos_theta, double freq);
	/// returns the complex gain for ideal microphone given alpha and beta
	static std::complex<double> microphoneDirectivity(double freq, double cos_theta, double alpha, double beta);
	/// returns the complex gain for ideal gradient microphone.
	static std::complex<double> gradientMicrophoneDirectivity(double freq, double cos_theta);


	/// coordinates of the microphone.
	double x;
	double y;
	double z;
	/// direction of the microphone.
	double direction;
	/// elevation of the microphone.
	double elevation;
};

#endif /* MICROPHONE_H_ */
