#ifndef GSCBEAMFORMER_H_
#define GSCBEAMFORMER_H_

#include "KinectConfig.h"
#include "SoundSourceLocalizer.h"

namespace Beam{
#define BM_N 16
#define BM_P 5
#define MC_Q 10
#define MC_L 16
	class GSCBeamformer {
	public:
		GSCBeamformer();
		~GSCBeamformer();
		void compute(float output[FRAME_SIZE], float input[][FRAME_SIZE], float angle, bool voice, float ref[FRAME_SIZE] = NULL);
	private:
		float m_input_prev[MAX_MICROPHONES][FRAME_SIZE];
		float m_d_prev[FRAME_SIZE];
		float m_y_prev[MAX_MICROPHONES][FRAME_SIZE];
		float m_bm[MAX_MICROPHONES][BM_N];
		float m_mc[MAX_MICROPHONES][MC_L];
	};
}

#endif /* GSCBEAMFORMER_H_ */
