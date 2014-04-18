SET(libbeam_sources
	Beamformer.h Beamformer.cpp
	Calibrator.h Calibrator.cpp
	Coords.h
	DSPFilter.h DSPFilter.cpp
	FFT.h FFT.cpp
	fftw3.h
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
	WavWriter.h WavWriter.cpp
)
ADD_LIBRARY(libbeam ${libbeam_sources})

TARGET_LINK_LIBRARIES(libbeam)

