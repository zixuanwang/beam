#include "beam\lib\Beamformer.h"
#include "beam\lib\FFT.h"
#include "beam\lib\Pipeline.h"
#include "beam\lib\Utils.h"
#include <iostream>

int main(int argc, char* argv[]){
	//// example to run the beamformer.
	//// sampling rate is 16000.
	Beam::WavReader wr("c:/dropbox/microsoft/beam/1.wav");
	char buf[50000];
	short* ptr = (short*)buf;
	int buf_filled;
	wr.read(buf, 50000, &buf_filled);
	for (int frame = 0; frame < 32; ++frame){
		std::stringstream ss;
		ss << "c:/users/danwa/desktop/" << frame << ".txt";
		std::ofstream out_stream(ss.str());
		std::vector<float> input[4]; // four channels in the time domain. each channel must have 512 samples.
		// fill up the input.
		for (int channel = 0; channel < 4; ++channel){
			for (int i = 0; i < 512; ++i){
				input[channel].push_back((float)(ptr[4 * i + channel]));
			}
		}
		std::vector<float> output(512, 0.f); // the output has 256 samples delay.
		std::vector<std::complex<float> > frequency_input[4]; // four channels in the frquency domain
		for (int channel = 0; channel < 4; ++channel){
			frequency_input[channel].assign(256, std::complex<float>(0.f, 0.f));
		}
		std::vector<std::complex<float> > beamformer_output(256); // 256 frequency bins
		Beam::FFT fft;
		for (int channel = 0; channel < 4; ++channel){
			fft.analyze(input[channel], frequency_input[channel]); // do FFT
		}
		// output to txt
		for (int channel = 0; channel < 4; ++channel){
			for (int bin = 0; bin < 256; ++bin){
				out_stream << frequency_input[channel][bin] << std::endl;
			}
		}
		Beam::Pipeline::instance()->preprocess(frequency_input); // phase compensation
		float angle;
		Beam::Pipeline::instance()->source_localize(frequency_input, &angle); // sound source localization & noise suppression
		Beam::Pipeline::instance()->smart_calibration(); // calibration
		Beam::Pipeline::instance()->beamforming(frequency_input, beamformer_output); // beamforming
		//Beam::Pipeline::instance()->postprocessing(beamformer_output); // frequency shifting
		fft.synthesize(beamformer_output, output); // do IFFT
		ptr += 256 * 4;
	}
	return 0;
}