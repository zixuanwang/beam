SET(libbeam_sources
	_kiss_fft_guts.h
	kiss_fft.h kiss_fft.c
	Math.h Math.cpp
	Microphone.h Microphone.cpp
	MicrophoneArray.h MicrophoneArray.cpp
	SoundSourceLocalizer.h SoundSourceLocalizer.cpp
	
)
ADD_LIBRARY(libbeam ${libbeam_sources})

