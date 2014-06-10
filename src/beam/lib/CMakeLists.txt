SET(libbeam_sources
	Beamformer.h Beamformer.cpp
	Calibrator.h Calibrator.cpp
	Coords.h
	DelaySumBeamformer.h DelaySumBeamformer.cpp
	DeReverb.h DeReverb.cpp
	DSPFilter.h DSPFilter.cpp
	FFT.h FFT.cpp
	GlobalConfig.h
	KinectConfig.h KinectConfig.cpp
	MicArrayDescriptor.h
	MicArrayWeights.h
	Microphone.h Microphone.cpp
	MVDRBeamformer.h MVDRBeamformer.cpp
	NoiseSuppressor.h NoiseSuppressor.cpp
	Pipeline.h Pipeline.cpp
	SoundSourceLocalizer.h SoundSourceLocalizer.cpp
	Tracker.h Tracker.cpp
	Utils.h
	WavReader.h WavReader.cpp
	WavWriter.h WavWriter.cpp
)
ADD_LIBRARY(libbeam ${libbeam_sources})

TARGET_LINK_LIBRARIES(libbeam) 

