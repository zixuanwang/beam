SET(libbeam_sources
	Coords.h
	DSPFilter.h DSPFilter.cpp
	fftw3.h
	GlobalConfig.h
	KinectConfig.h
	MCLT.h MCLT.cpp
	MicArrayDescriptor.h
	MicArrayWeights.h
	Microphone.h Microphone.cpp
	NoiseSuppressor.h NoiseSuppressor.cpp
	Pipeline.h Pipeline.cpp
	SoundSourceLocalizer.h SoundSourceLocalizer.cpp
	Tracker.h Tracker.cpp
	Utils.h
)
ADD_LIBRARY(libbeam ${libbeam_sources})

TARGET_LINK_LIBRARIES(libbeam)

