#include "WavWriter.h"

namespace Beam{
	WavWriter::WavWriter(const std::string& file, int sample_rate, short channels, short bit_per_sample, short format) : m_sample_rate(sample_rate), m_channels(channels), m_bit_per_sample(bit_per_sample), m_size(0){
		m_out_stream.open(file, std::ios::binary);
		m_out_stream.write("RIFF", 4);
		write<int>(m_out_stream, 36);
		m_out_stream.write("WAVE", 4);
		m_out_stream.write("fmt ", 4);
		write<int>(m_out_stream, 16);
		write<short>(m_out_stream, format);													// Format
		write<short>(m_out_stream, m_channels);											// Channels
		write<int>(m_out_stream, m_sample_rate);										// Sample Rate
		write<int>(m_out_stream, m_sample_rate * m_channels * m_bit_per_sample / 8);	// Byterate
		write<short>(m_out_stream, m_channels * m_bit_per_sample / 8);					// Frame size
		write<short>(m_out_stream, m_bit_per_sample);									// Bits per sample
		m_out_stream.write("data", 4);
	}
	
	WavWriter::~WavWriter(){
		m_out_stream.close();
	}

	void WavWriter::write(char* buf, int buf_size){
		// modify the header
		m_size += buf_size;
		m_out_stream.seekp(4, m_out_stream.beg);
		write<int>(m_out_stream, 36 + m_size);
		m_out_stream.seekp(40, m_out_stream.beg);
		m_out_stream.write((const char*)&m_size, sizeof(m_size));
		// modify the data
		m_out_stream.seekp(0, m_out_stream.end);
		m_out_stream.write((const char*)buf, buf_size);
	}
}
