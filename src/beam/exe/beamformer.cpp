#include "beam\lib\Beamformer.h"
#include "beam\lib\FFT.h"
#include "beam\lib\Pipeline.h"
#include "beam\lib\Utils.h"
#include <iostream>

std::string input_file;
std::string output_file;

void exit_with_help() {
	std::cout
		<< "Usage: beamformer input_file output_file\n";
	exit(1);
}

void parse_command_line(int argc, char* argv[]) {
	if (argc != 3)
		exit_with_help();
	input_file = argv[1];
	output_file = argv[2];
}

int main(int argc, char* argv[]){
	parse_command_line(argc, argv);
	Beam::WavReader reader(input_file);
	Beam::WavWriter writer(output_file, 16000, 1);
	Beam::FFT fft;
	int buf_size = FRAME_SIZE * 8;
	int output_buf_size = FRAME_SIZE * 2;
	char* buf = new char[buf_size];
	char* output_buf = new char[output_buf_size];
	short* output_ptr = (short*)output_buf;
	std::vector<float> input[MAX_MICROPHONES];
	std::vector<float> input_current[MAX_MICROPHONES];
	std::vector<float> input_prev[MAX_MICROPHONES];
	// initialize buffers
	for (int channel = 0; channel < MAX_MICROPHONES; ++channel){
		input[channel].assign(2 * FRAME_SIZE, 0.f);
		input_current[channel].assign(FRAME_SIZE, 0.f);
		input_prev[channel].assign(FRAME_SIZE, 0.f);
	}
	while (true){
		int buf_filled = 0;
		reader.read(buf, buf_size, &buf_filled);
		if (buf_filled < buf_size)
			break;
		reader.convert_format(input_current, buf, buf_size);
		for (int channel = 0; channel < MAX_MICROPHONES; ++channel){
			for (int i = 0; i < FRAME_SIZE; ++i){
				input[channel][i] = input_prev[channel][i];
			}
			for (int i = FRAME_SIZE; i < 2 * FRAME_SIZE; ++i){
				input[channel][i] = input_current[channel][i - FRAME_SIZE];
			}
		}
		std::vector<std::complex<float> > beamformer_output(FRAME_SIZE); // 256 frequency bins
		std::vector<float> output(2 * FRAME_SIZE, 0.f); // the output has 256 samples delay.
		std::vector<std::complex<float> > frequency_input[MAX_MICROPHONES]; // four channels in the frquency domain
		for (int channel = 0; channel < MAX_MICROPHONES; ++channel){
			frequency_input[channel].assign(FRAME_SIZE, std::complex<float>(0.f, 0.f));
		}
		for (int channel = 0; channel < MAX_MICROPHONES; ++channel){
			fft.analyze(input[channel], frequency_input[channel]); // do FFT
		}
		Beam::Pipeline::instance()->preprocess(frequency_input); // phase compensation
		float angle;
		Beam::Pipeline::instance()->source_localize(frequency_input, &angle); // sound source localization & noise suppression
		Beam::Pipeline::instance()->smart_calibration(); // calibration
		Beam::Pipeline::instance()->beamforming(frequency_input, beamformer_output); // beamforming
		Beam::Pipeline::instance()->postprocessing(beamformer_output); // frequency shifting
		fft.synthesize(beamformer_output, output); // do IFFT
		for (int i = 0; i < FRAME_SIZE; ++i){
			output_ptr[i] = (short)output[i];
		}
		writer.write(output_buf, FRAME_SIZE * 2);
		// copy half frame.
		for (int channel = 0; channel < MAX_MICROPHONES; ++channel){
			input_prev[channel].assign(input_current[channel].begin(), input_current[channel].end());
		}
	}
	delete[] buf;
	delete[] output_buf;
	std::cout << "conversion done..." << std::endl;
	return 0;
}