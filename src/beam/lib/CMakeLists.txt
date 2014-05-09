SET(libbeam_sources
	Beamformer.h Beamformer.cpp
	Calibrator.h Calibrator.cpp
	Coords.h
	DSPFilter.h DSPFilter.cpp
	FFT.h FFT.cpp
	GlobalConfig.h
	KinectConfig.h KinectConfig.cpp
	MCLT.h MCLT.cpp
	MicArrayDescriptor.h
	MicArrayWeights.h
	Microphone.h Microphone.cpp
	NoiseSuppressor.h NoiseSuppressor.cpp
	Pipeline.h Pipeline.cpp
	SoundSourceLocalizer.h SoundSourceLocalizer.cpp
	Tracker.h Tracker.cpp
	Utils.h
	WavReader.h WavReader.cpp
	WavWriter.h WavWriter.cpp
)
ADD_LIBRARY(libbeam ${libbeam_sources})

TARGET_LINK_LIBRARIES(libbeam ${CMAKE_PREFIX_PATH}/libfftw3f-3.lib ${CMAKE_PREFIX_PATH}/libfftw3-3.lib) 

