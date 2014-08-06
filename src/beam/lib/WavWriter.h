#ifndef WAVWRITER_H_
#define WAVWRITER_H_

#include <fstream>
#include <string>

namespace Beam{
	class WavWriter {
	public:
		WavWriter(const std::string& file, int sample_rate, short channels, short bit_per_sample, short format = 1);
		~WavWriter();
		void write(char* buf, int buf_size);
	private:
		std::ofstream m_out_stream;
		int m_sample_rate;
		short m_channels;
		int m_bit_per_sample;
		int m_size;
		template <typename T>
		void write(std::ofstream& stream, const T& t) {
			stream.write((const char*)&t, sizeof(T));
		}
	};
}

#endif /* WAVWRITER_H_ */
