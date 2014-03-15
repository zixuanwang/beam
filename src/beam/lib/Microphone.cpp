#include "Microphone.h"

Microphone::Microphone(){

}

Microphone::~Microphone(){

}

std::complex<double> Microphone::micRatio(double cos_theta, double freq){
	return microphoneDirectivity(freq, cos_theta, 1.f, 0.f);
}

std::complex<double> Microphone::microphoneDirectivity(double freq, double cos_theta, double alpha, double beta){
	return std::complex<double>(alpha, 0.0) + (gradientMicrophoneDirectivity(freq, cos_theta) * beta);
}

std::complex<double> Microphone::gradientMicrophoneDirectivity(double freq, double cos_theta){
	//	Parameters for typical ideal gradient microphone
	double delta = 0.005; // 5 millimeters delay
	double tau = delta / SOUND_SPEED;
	double filt_tau = 0.01;
	double r = 10e6;
	double c = filt_tau / r;
	std::complex<double> z_c, z_r, gain, z_c_1000, ratio_1000, gain_1000, ratio;
	std::complex<double> one(1.0, 0.0);
	if (freq < 1e-10){
		z_c = std::complex<double>(0.0, 1e-10);
	}
	else{
		z_c = std::complex<double>(0.0, 1.0 / (2.0 * PI * freq * c));
	}
	z_r = std::complex<double>(r, 0.0);
	gain = z_c / (z_c + z_r);
	z_c_1000 = std::complex<double>(0.0, 1.0 / (2.0 * PI * 1000.0 * c));
	gain_1000 = z_c_1000 / (z_c_1000 + z_r);
	double im = 2000.0 * PI * tau;
	ratio_1000 = one - (one / std::complex<double>(cos(im), sin(im)));
	ratio_1000 *= gain_1000;
	im = 2.0 * PI * freq * tau * cos_theta;
	ratio = one - (one / std::complex<double>(cos(im), sin(im)));
	ratio *= gain;
	ratio /= ratio_1000;
	return ratio;
}