#include "beam\lib\Beamformer.h"
#include "beam\lib\FFT.h"
#include "beam\lib\Pipeline.h"
#include "beam\lib\Utils.h"
#include <iostream>

int main(int argc, char* argv[]){
	//// example to run the beamformer.
	//// sampling rate is 16000.
	//// you need to fill up four input channels using captured data.
	//std::vector<float> input[4]; // four channels in the time domain. each channel must have 512 samples.
	//std::vector<float> output(512, 0.f); // the output has 256 samples delay.
	//std::vector<std::complex<float> > frequency_input[4]; // four channels in the frquency domain
	//for (int channel = 0; channel < 4; ++channel){
	//	frequency_input[channel].assign(256, std::complex<float>(0.f, 0.f));
	//}
	//std::vector<std::complex<float> > beamformer_output(256); // 256 frequency bins
	//Beam::FFT fft;
	//for (int channel = 0; channel < 4; ++channel){
	//	fft.analyze(input[channel], frequency_input[channel]); // do FFT
	//}
	//Beam::Pipeline::instance()->preprocess(frequency_input); // phase compensation
	//float angle;
	//Beam::Pipeline::instance()->source_localize(frequency_input, &angle); // sound source localization & noise suppression
	//Beam::Pipeline::instance()->smart_calibration(angle, frequency_input); // calibration
	//Beam::Pipeline::instance()->beamforming(frequency_input, beamformer_output); // beamforming
	//Beam::Pipeline::instance()->postprocessing(beamformer_output); // frequency shifting
	//fft.synthesize(beamformer_output, output); // do IFFT
	Beam::WavReader wr("c:/dropbox/microsoft/beam/1.wav");
	char buf[50000];
	short* ptr = (short*)buf;
	int buf_filled;
	wr.read(buf, 50000, &buf_filled);
	std::cout << (float)ptr[400] / 32767.f << std::endl;
	std::cout << (float)ptr[401] / 32767.f << std::endl;
	std::cout << (float)ptr[402] / 32767.f << std::endl;
	std::cout << (float)ptr[403] / 32767.f << std::endl;
	return 0;
}