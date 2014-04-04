#ifndef MICROPHONE_H_
#define MICROPHONE_H_

#include "GlobalConfig.h"
#include "Utils.h"

namespace Beam{
	class Microphone {
	public:
		Microphone(int _id, float _x, float _y, float _z, int _type, float _direction, float _elevation);
		~Microphone();
		/// microphone id.
		int id;
		/// coordinates of the microphone.
		float x;
		float y;
		float z;
		/// microphone type.
		int type;
		/// direction of the microphone.
		float direction;
		/// elevation of the microphone.
		float elevation;
		/// returns the complex gain for cardioid microphones.
		static std::complex<float> micRatio(float cos_theta, float freq);
		/// returns the complex gain for ideal microphone given alpha and beta
		static std::complex<float> microphoneDirectivity(float freq, float cos_theta, float alpha, float beta);
		/// returns the complex gain for ideal gradient microphone.
		static std::complex<float> gradientMicrophoneDirectivity(float freq, float cos_theta);
	};
}

#endif /* MICROPHONE_H_ */
