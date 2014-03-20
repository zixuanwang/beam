SET(libbeam_sources
	fftw3.h
	DSPFilter.h DSPFilter.cpp
	GlobalConfig.h
	MCLT.h MCLT.cpp
	Microphone.h Microphone.cpp
	MicrophoneArray.h MicrophoneArray.cpp
	NoiseSuppressor.h NoiseSuppressor.cpp
	Pipeline.h Pipeline.cpp
	SoundSourceLocalizer.h SoundSourceLocalizer.cpp
	Tracker.h Tracker.cpp
	Utils.h
)
ADD_LIBRARY(libbeam ${libbeam_sources})

TARGET_LINK_LIBRARIES(libbeam)

