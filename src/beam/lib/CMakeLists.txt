SET(libbeam_sources
	fftw3.h
	Math.h Math.cpp
	MCLT.h MCLT.cpp
	Microphone.h Microphone.cpp
	MicrophoneArray.h MicrophoneArray.cpp
	SoundSourceLocalizer.h SoundSourceLocalizer.cpp
	
)
ADD_LIBRARY(libbeam ${libbeam_sources})

TARGET_LINK_LIBRARIES(libbeam)

