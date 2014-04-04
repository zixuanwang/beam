#include "Microphone.h"

namespace Beam{
	Microphone::Microphone(int _id, float _x, float _y, float _z, int _type, float _direction, float _elevation) : id(_id), x(_x), y(_y), z(_z), type(_type), direction(_direction), elevation(_elevation){
	}

	Microphone::~Microphone(){

	}

	std::complex<float> Microphone::micRatio(float cos_theta, float freq){
		return microphoneDirectivity(freq, cos_theta, 0.5f, 0.5f);
	}

	std::complex<float> Microphone::microphoneDirectivity(float freq, float cos_theta, float alpha, float beta){
		return std::complex<float>(alpha, 0.f) + (gradientMicrophoneDirectivity(freq, cos_theta) * beta);
	}

	std::complex<float> Microphone::gradientMicrophoneDirectivity(float freq, float cos_theta){
		//	Parameters for typical ideal gradient microphone
		float delta = 0.005f; // 5 millimeters delay
		float tau = (float)(delta / SOUND_SPEED);
		float filt_tau = 0.01f;
		float r = 10e6f;
		float c = filt_tau / r;
		std::complex<float> z_c, z_r, gain, z_c_1000, ratio_1000, gain_1000, ratio;
		std::complex<float> one(1.f, 0.f);
		if (freq < 1e-10f){
			z_c = std::complex<float>(0.f, 1e-10f);
		}
		else{
			z_c = std::complex<float>(0.f, 1.f / ((float)(2.f * PI * freq * c)));
		}
		z_r = std::complex<float>(r, 0.f);
		gain = z_c / (z_c + z_r);
		z_c_1000 = std::complex<double>(0.f, 1.f / (2.f * PI * 1000.f * c));
		gain_1000 = z_c_1000 / (z_c_1000 + z_r);
		float im = (float)(2000.f * PI * tau);
		ratio_1000 = one - (one / std::complex<float>(cosf(im), sinf(im)));
		ratio_1000 *= gain_1000;
		im = (float)(2.f * PI * freq * tau * cos_theta);
		ratio = one - (one / std::complex<float>(cosf(im), sinf(im)));
		ratio *= gain;
		ratio /= ratio_1000;
		return ratio;
	}
}