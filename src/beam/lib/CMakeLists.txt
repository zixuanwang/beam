SET(libbeam_sources
	Math.h Math.cpp
	Microphone.h Microphone.cpp
	MicrophoneArray.h MicrophoneArray.cpp
	SoundSourceLocalizer.h SoundSourceLocalizer.cpp
	
)
ADD_LIBRARY(libbeam ${libbeam_sources})

