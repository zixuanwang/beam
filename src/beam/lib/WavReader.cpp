#include "WavReader.h"
#include <iostream>

namespace Beam{
	WavReader::WavReader(const std::string& file){
		m_in_stream.open(file, std::ios::binary);
		if (m_in_stream.good()){
			m_in_stream.seekg(22, m_in_stream.beg);
			m_in_stream.read((char*)&m_channels, sizeof(m_channels));
			m_in_stream.read((char*)&m_sample_rate, sizeof(m_sample_rate));
			m_in_stream.seekg(34, m_in_stream.beg);
			m_in_stream.read((char*)&m_bit_per_sample, sizeof(m_bit_per_sample));
			m_in_stream.seekg(40, m_in_stream.beg);
			m_in_stream.read((char*)&m_size, sizeof(m_size));
			m_in_stream.seekg(0, m_in_stream.end);
			m_data_length = (long long)m_in_stream.tellg();
			m_data_length -= 44;
			m_in_stream.seekg(44, m_in_stream.beg);
		}
	}

	WavReader::~WavReader(){
		if (m_in_stream.good()){
			m_in_stream.close();
		}
	}

	void WavReader::read(char* buf, int buf_size, int* filled_size){
		if (m_in_stream.good()){
			if (m_data_length < (long long)buf_size){
				m_in_stream.read(buf, m_data_length);
				*filled_size = (int)m_data_length;
				m_data_length = 0;
			}
			else{
				m_in_stream.read(buf, buf_size);
				m_data_length -= buf_size;
				*filled_size = buf_size;
			}
		}
	}

	//void WavReader::convert_format(std::vector<float>* input, char* buf, int buf_size){
	//	short* ptr = (short*)(buf);
	//	int len = buf_size / m_bit_per_sample * 8 / m_channels;
	//	for (int channel = 0; channel < m_channels; ++channel){
	//		input[channel].assign(len, 0.f);
	//		for (int bin = 0; bin < len; ++bin){
	//			input[channel][bin] = (float)(ptr[m_channels * bin + channel]);
	//		}
	//	}
	//}

	void WavReader::convert_format(float input[][FRAME_SIZE], char* buf, int buf_size){
		if (m_bit_per_sample == 16){
			short* ptr = (short*)(buf);
			int len = buf_size / m_bit_per_sample * 8 / m_channels;
			for (int channel = 0; channel < m_channels; ++channel){
				for (int bin = 0; bin < len; ++bin){
					input[channel][bin] = (float)(ptr[m_channels * bin + channel]) / SHRT_MAX;
				}
			}
		}
		if (m_bit_per_sample == 32){
			// for new kinect
			int* ptr = (int*)(buf);
			int len = buf_size / m_bit_per_sample * 8 / m_channels;
			for (int channel = 0; channel < m_channels; ++channel){
				for (int bin = 0; bin < len; ++bin){
					// switch channels because of different geometry.
					input[MAX_MICROPHONES - 1 - channel][bin] = (float)ptr[m_channels * bin + channel] / INT_MAX;			
				}
			}
		}
	}

	int WavReader::get_channels(){
		return (int)m_channels;
	}

	int WavReader::get_bit_per_sample(){
		return (int)m_bit_per_sample;
	}

	int WavReader::swap_int32(int val){
		val = ((val << 8) & 0xFF00FF00) | ((val >> 8) & 0xFF00FF);
		return (val << 16) | ((val >> 16) & 0xFFFF);
	}
}
