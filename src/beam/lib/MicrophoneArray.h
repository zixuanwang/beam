#ifndef MICROPHONEARRAY_H_
#define MICROPHONEARRAY_H_

#include "Microphone.h"

#define MAX_MICROPHONES 4

class MicrophoneArray {
public:
	MicrophoneArray();
	~MicrophoneArray();
	Microphone mic_array[MAX_MICROPHONES];

};

#endif /* MICROPHONEARRAY_H_ */
