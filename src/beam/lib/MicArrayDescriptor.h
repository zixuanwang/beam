#ifndef MICARRAYDESCRIPTOR_H_
#define MICARRAYDESCRIPTOR_H_

#include "Microphone.h"

namespace Beam{
	class MicArrayDescriptor {
	public:
		char manifacturer[64]; // manufacturer string
		char model[64]; // model string
		int mic_array_type; // type of the microphone array
		// work volume definition:
		float work_vert_angle_beg; // vertical angle begin, rad
		float work_vert_angle_end; // vertical angle end, rad
		float work_hor_angle_beg; // horizontal angle beg, rad
		float work_hor_angle_end; // horizontal angle end, rad
		// Work frequency band definition:
		float freq_low; // Frequency range - lower, Hz
		float freq_high; // Frequency range - upper, Hz
		// Microphones descriptors:
		int num_mics; // number of microphones in the array
		Microphone mic[MAX_MICROPHONES];
	};
}



#endif /* MICARRAYDESCRIPTOR_H_ */
