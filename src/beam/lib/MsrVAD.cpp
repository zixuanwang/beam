#include "MsrVAD.h"

namespace Beam{
	MsrVAD::MsrVAD()
	{
		m_VAD_TAUN = 0.4f;
		m_VAD_TAUS = 2.681769f;
		m_VAD_SNR_SMOOTHER = 0.941634f;
		m_VAD_A01 = 0.05f;
		m_VAD_A10 = 0.99f;
		m_VAD_A01f = 0.01f;
		m_VAD_A10f = 0.870009f;
		m_MecBegBin = std::max(2, (int)(200 * (SAMPLE_RATE / 1000) / FRAME_SIZE / 2));
		m_MecEndBin = std::min(FRAME_SIZE - 1, (int)(7200 * (SAMPLE_RATE / 1000) / FRAME_SIZE / 2));
		m_MecFrameDuration = ((float)FRAME_SIZE / SAMPLE_RATE);

		m_LogLikelihoodBegBin = 0.045f * FRAME_SIZE;
		m_LogLikelihoodEndBin = 0.65f * FRAME_SIZE;

		std::fill(q30MasterSpeechPresenceProb, q30MasterSpeechPresenceProb + FRAME_SIZE, 0);
		std::fill(MasterSignalPowerOverNoiseModel, MasterSignalPowerOverNoiseModel + FRAME_SIZE, 0.f);
		std::fill(priorSNR, priorSNR + FRAME_SIZE, 0.f);
		std::fill(posteriorSNR_NS, posteriorSNR_NS + FRAME_SIZE, 0.f);
		std::fill(speechPresenceLR, speechPresenceLR + FRAME_SIZE, 0.f);
		std::fill(speechPresenceProb, speechPresenceProb + FRAME_SIZE, 0.f);
		std::fill(SpeechModel, SpeechModel + FRAME_SIZE, 1.f);
		std::fill(NoiseModel, NoiseModel + FRAME_SIZE, 0.f);
		std::fill(SignalSpecIn, SignalSpecIn + 2 * FRAME_SIZE, 0.f);
		std::fill(SignalPower, SignalPower + FRAME_SIZE, 0.f);

		FrameLR = 0.0f;
		m_fFramePresProb = 0;
		nFrame = 1;
		m_fEnergy = 0.0f;
		m_fSNR = 1.0f;
	}

	MsrVAD::~MsrVAD()
	{
	}

	void MsrVAD::process(float* fft_ptr){
		const float a00 = 1.0f - m_VAD_A01;
		const float a11 = 1.0f - m_VAD_A10;
		const float a00f = 1.0f - m_VAD_A01f;
		const float a11f = 1.0f - m_VAD_A10f;

		//This XBox code expects the data in the [-1.0,1.0] range so scale appropiately
		//const FLOAT32 toFloat = float(1u << (ctrVAD ? 25 : 23));
		//const float toFloat = float(1u << 15);

		//for (unsigned int u = 0; u < TWO_FRAME_SIZE; ++u)
		//{
		//	SignalSpecIn[u] = (float)fft_ptr[u] / toFloat;
		//}

		// per bin and per frame soft VAD
		float likemean = 0.0f;
		float likelogmean = 0.0f;
		for (unsigned int u = m_MecBegBin; u <= m_MecEndBin; u++) { // optimized for non-zero frequency bins
			// calculate signal power
			float MicRe = fft_ptr[u];
			float MicIm = fft_ptr[TWO_FRAME_SIZE - u];
			SignalPower[u] = MicRe * MicRe + MicIm * MicIm;

			// update the prior and posterios SNRs
			posteriorSNR_NS[u] = SignalPower[u] / (NoiseModel[u] + DIV_BY_ZERO_PREVENTION);
			float MLPriorSNR = SpeechModel[u] / (NoiseModel[u] + DIV_BY_ZERO_PREVENTION) - 1.0f;
			if (MLPriorSNR < 0.0f) MLPriorSNR = 0.0f;
			priorSNR[u] = m_VAD_SNR_SMOOTHER * priorSNR[u] + (1.0f - m_VAD_SNR_SMOOTHER) * MLPriorSNR;

			// Speech presence likelihood ratio
			float vRatio = posteriorSNR_NS[u] * priorSNR[u] / (priorSNR[u] + 1.0f);
			if (vRatio > 700.0f) vRatio = 700.0f;
			float likelihoodratio = 1.0f / (priorSNR[u] + 1.0f) * expf(vRatio);
			if (likelihoodratio > 1000.0f) likelihoodratio = 1000.0f;
			else if (likelihoodratio < LOG_LIKELIHOOD_MINVAL) likelihoodratio = LOG_LIKELIHOOD_MINVAL;

			// Smooth speech presence probability per bin
			speechPresenceLR[u] = likelihoodratio * (m_VAD_A01 + a11 * speechPresenceLR[u]) / (a00 + m_VAD_A10 * speechPresenceLR[u]); // HMM for changing the state
			speechPresenceProb[u] = speechPresenceLR[u] / (1.0f + speechPresenceLR[u]);
			if (speechPresenceProb[u] < 0.0f) speechPresenceProb[u] = 0.0f;
			else if (speechPresenceProb[u] > 1.0f) speechPresenceProb[u] = 1.0f;

			// Note that likelogmean is only calculated on a subset of the frequency bins
			if (u >= m_LogLikelihoodBegBin && u <= m_LogLikelihoodEndBin) {
				likemean += likelihoodratio;
				likelogmean += logf(likelihoodratio);
			}
		}
		likemean /= (m_LogLikelihoodEndBin - m_LogLikelihoodBegBin + 1);
		likelogmean /= (m_LogLikelihoodEndBin - m_LogLikelihoodBegBin + 1);

		// Speech presence likelihood ratio for the frame
		float curFrameLR = expf(likelogmean);
		curFrameLR = 0.8f * curFrameLR + (1 - 0.8f) * likemean;
		if (curFrameLR > 1000.0f) curFrameLR = 1000.0f;

		// Smooth speech presence probability per frame
		FrameLR = curFrameLR * (m_VAD_A01f + a11f * FrameLR) / (a00f + m_VAD_A10f * FrameLR); // HMM for changing the state
		float framePresProb = FrameLR / (1.0f + FrameLR);
		if (framePresProb < 0.0f) framePresProb = 0.0f;
		else if (framePresProb > 1.0f) framePresProb = 1.0f;
		m_fFramePresProb = framePresProb;

		// precise noise model
		if (nFrame > 1) {
			if (nFrame < m_VAD_TAUN / m_MecFrameDuration) {
				for (unsigned int u = m_MecBegBin; u <= m_MecEndBin; u++) {
					float alphaN = 1.0f / nFrame;
					NoiseModel[u] = (1.0f - alphaN) * NoiseModel[u] + alphaN * SignalPower[u];
				}
			}
			else {
				for (unsigned int u = m_MecBegBin; u <= m_MecEndBin; u++) {
					float alphaN = (1.0f - speechPresenceProb[u]) * (1.0f - framePresProb) * m_MecFrameDuration / m_VAD_TAUN;
					NoiseModel[u] = (1.0f - alphaN) * NoiseModel[u] + alphaN * SignalPower[u];
				}
			}

			// update the speech model
			for (unsigned int u = m_MecBegBin; u <= m_MecEndBin; u++) {
				float alphaS = speechPresenceProb[u] * framePresProb * m_MecFrameDuration / m_VAD_TAUS;
				SpeechModel[u] = (1.0f - alphaS) * SpeechModel[u] + alphaS * SignalPower[u];
			}
		}
		nFrame++;

		// Update the prior SNR
		m_fEnergy = 0.0f;
		float fSNRam = 0.0f;
		float fSNR = 0.0f;
		for (unsigned int u = m_MecBegBin; u <= m_MecEndBin; u++) {
			if (u >= m_LogLikelihoodBegBin && u <= m_LogLikelihoodEndBin) {
				fSNR = SignalPower[u] / (NoiseModel[u] + 1e-10f);
				fSNR = std::max(1.0f, std::min(20000.0f, fSNR));
				fSNRam += fSNR;
			}

			// Compute the suppression rule - simple Wiener
			float Gain = priorSNR[u] / (priorSNR[u] + 1.0f);
			// update the prior SNR
			MasterSignalPowerOverNoiseModel[u] = SignalPower[u] / (NoiseModel[u] + DIV_BY_ZERO_PREVENTION);
			float alphaS = speechPresenceProb[u] * framePresProb * m_MecFrameDuration / m_VAD_TAUS;
			priorSNR[u] = (1.0f - alphaS) * priorSNR[u] + alphaS *  Gain * Gain * MasterSignalPowerOverNoiseModel[u];
			m_fEnergy += SignalPower[u];
			q30MasterSpeechPresenceProb[u] = (int)((float)(1u << 30) * framePresProb * speechPresenceProb[u]);
		}
		m_fEnergy /= (m_MecEndBin - m_MecBegBin + 1);
		fSNRam /= (m_MecEndBin - m_MecBegBin + 1);
		float beta = m_MecFrameDuration*framePresProb / 10.0f;  // time constant in seconds
		m_fSNR = (1.0f - beta) * m_fSNR + beta * fSNRam;
	}
}
