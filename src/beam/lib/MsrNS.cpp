#include "MsrNS.h"

namespace Beam{
	MsrNS::MsrNS()
	{
		std::fill(piPriorSNR, piPriorSNR + FRAME_SIZE, 1.f);
		std::fill(piGain, piGain + FRAME_SIZE, 1.f);
		std::fill(piPosteriorSNR_NS, piPosteriorSNR_NS + FRAME_SIZE, 0.f);
		std::fill(piMLPriorSNR_NS, piMLPriorSNR_NS + FRAME_SIZE, 0.f);
		std::fill(piPriorSNR_NS, piPriorSNR_NS + FRAME_SIZE, 0.f);
	}

	MsrNS::~MsrNS()
	{
	}

	void MsrNS::process(float* fft_ptr, MsrVAD* pVAD, bool enableNS){
		int* pMasterSpeechPresenceProb = pVAD->GetMasterSpeechPresenceProbQ30();
		float* pMasterSignalPowerOverNoiseModel = pVAD->GetMasterSignalPowerOverNoiseModel();

		unsigned int mecBegBin = pVAD->GetBeginBin();
		unsigned int mecEndBin = pVAD->GetEndBin();

		// Floating point
		for (unsigned int i = mecBegBin; i <= mecEndBin; i++)
		{
			// Compute prior and posterior SNRs
			piPosteriorSNR_NS[i] = pMasterSignalPowerOverNoiseModel[i];
			piMLPriorSNR_NS[i] = piPosteriorSNR_NS[i] - 1.0f;
			if (piMLPriorSNR_NS[i] < 0.0f)
				piMLPriorSNR_NS[i] = 0.0f;

			piPriorSNR_NS[i] = (NS_SNR_SMOOTHER_NS_F * piPriorSNR[i]) + (ONEMINUS_SNR_SMOOTHER_NS_F * piMLPriorSNR_NS[i]);

			// Compute the suppression rule
			float fHGain = piPriorSNR_NS[i];

			if (fHGain < MINSNR_F)
			{
				fHGain = MINSNR_F;
			}
			else if (fHGain > MAXSNR_F)
			{
				fHGain = MAXSNR_F;
			}

			// Hgain
			fHGain = fHGain / (1.0f + fHGain);

			// Limit to -60dB
			if (fHGain < .001f)
			{
				fHGain = .001f;
			}

			float fHGainClipped = fHGain;
			if (fHGainClipped > 1.0f)
				fHGainClipped = 1.0f;

			fHGain = NS_MINGAIN_F + (fHGainClipped * ONEMINUS_MINGAIN_F);

			// Compute UncertainGain       
			float fUncertainGain = ((ONEMINUS_MINGAINU_F / float(1u << 30)) * pMasterSpeechPresenceProb[i]) + NS_MINGAINU_F;
			fUncertainGain *= fHGain;

			// Smooth the gain - calculate final NS gain
			float fGain = (ONEMINUS_BETA_F * piGain[i]) + (NS_BETA_F * fUncertainGain);  //(1.0-beta) * Gain[i] + beta * UncertainGain

			// Apply the gain to pSpec
			if (enableNS)
			{
				fft_ptr[i] *= fGain;                 // Re
				fft_ptr[TWO_FRAME_SIZE - i] *= fGain;
			}

			// Update gain for next frame
			piGain[i] = fGain;

			// Update prior SNR for the next frame
			piPriorSNR[i] = fGain * fGain * pMasterSignalPowerOverNoiseModel[i];
		}
	}
}
