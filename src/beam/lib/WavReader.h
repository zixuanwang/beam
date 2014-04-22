#ifndef WAVREADER_H_
#define WAVREADER_H_

#include <fstream>
#include <string>

namespace Beam{
	class WavReader {
	public:
		WavReader(const std::string& file);
		~WavReader();
		void read(char* buf, int buf_size, int* filled_size);
	private:
		std::ifstream m_in_stream;
		int m_sample_rate;
		short m_channels;
		int m_size;
	};
}

#endif /* WAVREADER_H_ */
