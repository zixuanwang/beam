#ifndef WAVREADER_H_
#define WAVREADER_H_

#include <fstream>
#include <string>
#include <vector>
#include "GlobalConfig.h"

namespace Beam{
	class WavReader {
	public:
		WavReader(const std::string& file);
		~WavReader();
		void read(char* buf, int buf_size, int* filled_size);
		/// convert the input buffer to appropriate format.
		//void convert_format(std::vector<float>* input, char* buf, int buf_size);
		void convert_format(float input[][FRAME_SIZE], char* buf, int buf_size);
		int get_channels();
		int get_bit_per_sample();
		int swap_int32(int val);
	private:
		std::ifstream m_in_stream;
		int m_sample_rate;
		short m_channels;
		short m_bit_per_sample;
		int m_size;
		long long m_data_length;
	};
}

#endif /* WAVREADER_H_ */
