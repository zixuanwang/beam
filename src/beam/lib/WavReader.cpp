#include "WavReader.h"
#include <iostream>

namespace Beam{
	WavReader::WavReader(const std::string& file){
		m_in_stream.open(file, std::ios::binary);
		if (m_in_stream.good()){
			m_in_stream.seekg(22, m_in_stream.beg);
			m_in_stream.read((char*)&m_channels, sizeof(m_channels));
			m_in_stream.read((char*)&m_sample_rate, sizeof(m_sample_rate));
			m_in_stream.seekg(40, m_in_stream.beg);
			m_in_stream.read((char*)&m_size, sizeof(m_size));
			m_in_stream.seekg(0, m_in_stream.end);
			int length = (int)m_in_stream.tellg();
			length -= 44;
			m_size = m_size < length ? m_size : length;
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
			if (m_size < buf_size){
				m_in_stream.read(buf, m_size);
				m_size = 0;
				*filled_size = m_size;
			}
			else{
				m_in_stream.read(buf, buf_size);
				m_size -= buf_size;
				*filled_size = buf_size;
			}
		}
	}

	void WavReader::convert_format(std::vector<float>* input, char* buf, int buf_size){
		short* ptr = (short*)(buf);
		int len = buf_size / 2 / MAX_MICROPHONES;
		for (int channel = 0; channel < MAX_MICROPHONES; ++channel){
			input[channel].assign(len, 0.f);
			for (int bin = 0; bin < len; ++bin){
				input[channel][bin] = (float)(ptr[MAX_MICROPHONES * bin + channel]);
			}
		}
	}
}
