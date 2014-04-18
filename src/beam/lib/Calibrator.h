#ifndef CALIBRATOR_H_
#define CALIBRATOR_H_

#include <cfloat>
#include "DSPFilter.h"
#include "KinectConfig.h"

namespace Beam{
	class Calibrator {
	public:
		Calibrator();
		~Calibrator();
		float calibrate(float sound_source, std::vector<std::complex<float> >* input, std::complex<float> persistent_gains[MAX_MICROPHONES][MAX_GAIN_SUBBANDS]);
	private:
		float m_coordinates[MAX_MICROPHONES];
		float m_coeff[2];
		std::vector<std::complex<float> > m_frequency_filter[MAX_GAIN_SUBBANDS];
		std::vector<std::complex<float> > m_working_frequency;
	};
}

#endif /* CALIBRATOR_H_ */
