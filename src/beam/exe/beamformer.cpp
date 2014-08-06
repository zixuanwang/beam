#include "beam/lib/Beamformer.h"
#include "beam/lib/FFT.h"
#include "beam/lib/MCLT.h"
#include "beam/lib/Pipeline.h"
#include "beam/lib/Utils.h"
#include <iostream>

std::string input_file;
std::string output_file;

void exit_with_help() {
	std::cout << "Usage: beamformer input_file output_file\n";
	exit(1);
}

void parse_command_line(int argc, char* argv[]) {
	if (argc != 3)
		exit_with_help();
	input_file = argv[1];
	output_file = argv[2];
}

int main(int argc, char* argv[]) {
	parse_command_line(argc, argv);
	Beam::WavReader reader(input_file);
	int channels = reader.get_channels();
	Beam::WavWriter writer(output_file, 16000, 1, 16);
	int buf_size = FRAME_SIZE * channels * 2;
	int output_buf_size = FRAME_SIZE * 2;
	char* buf = new char[buf_size];
	char* output_buf = new char[output_buf_size];
	short* output_ptr = (short*) output_buf;
	float input[MAX_MICROPHONES][FRAME_SIZE] = { 0.f };
	float output[FRAME_SIZE] = { 0.f };
	while (true) {
		int buf_filled = 0;
		reader.read(buf, buf_size, &buf_filled);
		if (buf_filled < buf_size)
			break;
		reader.convert_format(input, buf, buf_size);
		// this is the key step in the beamformer.
		// input are 4 channels. each channel contains 256 float numbers.
		// output is 1 channel. it contains 256 float numbers.
		Beam::Pipeline::instance()->process(input, output); 
		for (int i = 0; i < FRAME_SIZE; ++i) {
			output_ptr[i] = (short) output[i];
		}
		writer.write(output_buf, FRAME_SIZE * 2);
	}
	delete[] buf;
	delete[] output_buf;
	std::cout << "conversion done..." << std::endl;
	return 0;
}
