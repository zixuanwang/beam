#include "WavReader.h"

namespace Beam{
	WavReader::WavReader(const std::string& file){
		m_in_stream.open(file, std::ios::binary);
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

	WavReader::~WavReader(){
		if (m_in_stream.good()){
			m_in_stream.close();
		}
	}

	void WavReader::read(char* buf, int buf_size, int* filled_size){
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
