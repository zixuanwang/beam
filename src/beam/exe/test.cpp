#include "beam\lib\Beamformer.h"
#include "beam\lib\FFT.h"
#include "beam\lib\MCLT.h"
#include "beam\lib\Pipeline.h"
#include "beam\lib\Utils.h"
#include <chrono>
#include <iostream>

int main(int argc, char* argv[]){
	std::chrono::high_resolution_clock::time_point p1, p2;
	p1 = std::chrono::high_resolution_clock::now();
	std::complex<float> a(1.f, -3.f);
	std::complex<float> b(2.f, 2.3f);
	//std::complex<float> c;
	float sum = 0.f;
	for (int i = 0; i < 100000000; ++i){
		//float re = a.real();
		//float im = a.imag();
		//sum += sqrtf(re * re + im * im);
		//std::cout << re << std::endl;
		//a += b;
		//sum += std::norm(a);
		//std::complex<float> d = a;
		//d += b;
		//sum += std::norm(a+b);
		std::complex<float> c = a + b;
		sum += std::norm(c);
	}
	//float f1 = FLT_MAX;
	//std::cout << f1 << std::endl;
	//float f2 = f1 * f1;
	//float f3 = sqrtf(f2);
	//std::cout << f2 << std::endl;
	p2 = std::chrono::high_resolution_clock::now();
	std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(p2 - p1).count() << std::endl;
	std::cout << sum << std::endl;
	// example to run the beamformer.
	// sampling rate is 16000.
	//Beam::WavReader wr("c:/dropbox/microsoft/beam/1.wav");
	//int len = 4000000;
	//char* buf = new char[len];
	//short* ptr = (short*)buf;
	//int buf_filled;
	//wr.read(buf, len, &buf_filled);
	//for (int frame = 0; frame < 1280; ++frame){
	//	std::stringstream ss;
	//	ss << "c:/users/danwa/desktop/" << frame << ".txt";
	//	std::ofstream out_stream(ss.str());
	//	std::vector<float> input[4]; // four channels in the time domain. each channel must have 512 samples.
	//	// fill up the input.
	//	for (int channel = 0; channel < 4; ++channel){
	//		input[channel].assign(512, 0.f);
	//		for (int bin = 0; bin < 512; ++bin){
	//			input[channel][bin] = (float)(ptr[4 * bin + channel]);
	//		}
	//	}
	//	std::vector<float> output(512, 0.f); // the output has 256 samples delay.
	//	std::vector<std::complex<float> > frequency_input[4]; // four channels in the frquency domain
	//	for (int channel = 0; channel < 4; ++channel){
	//		frequency_input[channel].assign(256, std::complex<float>(0.f, 0.f));
	//	}
	//	std::vector<std::complex<float> > beamformer_output(256); // 256 frequency bins
	//	Beam::FFT fft;
	//	for (int channel = 0; channel < 4; ++channel){
	//		fft.analyze(input[channel], frequency_input[channel]); // do FFT
	//	}
	//	// output to txt
	//	for (int channel = 0; channel < 4; ++channel){
	//		for (int bin = 0; bin < 256; ++bin){
	//			out_stream << frequency_input[channel][bin] << std::endl;
	//		}
	//	}
	//	//Beam::Pipeline::instance()->preprocess(frequency_input); // phase compensation
	//	float angle;
	//	Beam::Pipeline::instance()->source_localize(frequency_input, &angle); // sound source localization & noise suppression
	//	//Beam::Pipeline::instance()->smart_calibration(); // calibration
	//	Beam::Pipeline::instance()->beamforming(frequency_input, beamformer_output); // beamforming
	//	for (auto& v : beamformer_output){
	//		std::cout << v << std::endl;
	//	}
	//	std::cout << frame << std::endl;
	//	//Beam::Pipeline::instance()->postprocessing(beamformer_output); // frequency shifting
	//	fft.synthesize(beamformer_output, output); // do IFFT
	//	ptr += 256 * 4;
	//}
	//delete[] buf;
	return 0;
}