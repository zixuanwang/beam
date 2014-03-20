#include "MicrophoneArray.h"

namespace Beam{
	MicrophoneArray::MicrophoneArray(){
		// default is the kinect microphone
		mic_array[0].x = 0.113f;
		mic_array[0].y = 0.02f;
		mic_array[0].z = 0.f;
		mic_array[0].direction = 0.f;
		mic_array[0].elevation = 0.f;
		mic_array[1].x = -0.036f;
		mic_array[1].y = 0.02f;
		mic_array[1].z = 0.f;
		mic_array[1].direction = 0.f;
		mic_array[1].elevation = 0.f;
		mic_array[2].x = -0.076f;
		mic_array[2].y = 0.02f;
		mic_array[2].z = 0.f;
		mic_array[2].direction = 0.f;
		mic_array[2].elevation = 0.f;
		mic_array[3].x = -0.113f;
		mic_array[3].y = 0.02f;
		mic_array[3].z = 0.f;
		mic_array[3].direction = 0.f;
		mic_array[3].elevation = 0.f;
	}

	MicrophoneArray::~MicrophoneArray(){

	}
}
