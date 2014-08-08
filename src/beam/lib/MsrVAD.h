#ifndef MSRVAD_H_
#define MSRVAD_H_

#include <algorithm>
#include "GlobalConfig.h"
#include "Utils.h"

namespace Beam{
#define DIV_BY_ZERO_PREVENTION 1e-10f    // Very small number used in y = 1/(x) type calculations to prevent dividing by zero.
#define LOG_LIKELIHOOD_MINVAL  1e-10f    // minimum allowed value for likelihoodratio
	class MsrVAD
	{
	public:
		MsrVAD();
		~MsrVAD();
		void process(float* fft_ptr);
		int* GetMasterSpeechPresenceProbQ30() { return &q30MasterSpeechPresenceProb[0]; }
		float* GetMasterSignalPowerOverNoiseModel() { return &MasterSignalPowerOverNoiseModel[0]; }
		float GetSNR() { return m_fSNR; }
		unsigned int GetBeginBin() { return m_MecBegBin; }
		unsigned int GetEndBin() { return m_MecEndBin; }
	private:
		float m_VAD_TAUN;
		float m_VAD_TAUS;
		float m_VAD_SNR_SMOOTHER;
		float m_VAD_A01;
		float m_VAD_A10;
		float m_VAD_A01f;
		float m_VAD_A10f;
		unsigned int m_MecBegBin;
		unsigned int m_MecEndBin;
		float  m_MecFrameDuration;
		float m_LogLikelihoodBegBin;
		float m_LogLikelihoodEndBin;

		int q30MasterSpeechPresenceProb[FRAME_SIZE];
		float MasterSignalPowerOverNoiseModel[FRAME_SIZE];
		float priorSNR[FRAME_SIZE];
		float posteriorSNR_NS[FRAME_SIZE];
		float speechPresenceLR[FRAME_SIZE]; // per bin soft VAD
		float speechPresenceProb[FRAME_SIZE];
		float SpeechModel[FRAME_SIZE];
		float NoiseModel[FRAME_SIZE];
		float SignalSpecIn[2 * FRAME_SIZE];
		float SignalPower[FRAME_SIZE];

		float FrameLR;
		float m_fFramePresProb;
		float m_fEnergy;
		float m_fSNR;
		unsigned int nFrame;

	};
}

#endif /* MSRVAD_H_ */
