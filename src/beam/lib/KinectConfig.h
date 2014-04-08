#ifndef KINECTCONFIG_H_
#define KINECTCONFIG_H_

#include "MicArrayDescriptor.h"
#include "MicArrayWeights.h"

namespace Beam{
	class KinectConfig{
	public:
		KinectConfig();
		~KinectConfig();
		static MicArrayWeights kinect_weights;
		static MicArrayDescriptor kinect_descriptor;
		static float frequency_bands[MAX_GAIN_SUBBANDS][3];
		static RCoords kinect_beams[11];
		static float kinect_frequencies[256];
		static std::complex<float> kinect_comp_weights[11264];
		static std::complex<float> kinect_dd[11][256][4];
	};
}

#endif /* KINECTCONFIG_H_ */
