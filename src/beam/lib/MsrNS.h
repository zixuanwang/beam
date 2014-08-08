#ifndef MSRNS_H_
#define MSRNS_H_

#include "GlobalConfig.h"
#include "MsrVAD.h"
#include "Utils.h"

namespace Beam{
#define NS_SNR_SMOOTHER_NS_F 0.79f
#define NS_MINGAIN_F 0.4f
#define NS_MINGAINU_F 0.44f
#define NS_BETA_F 0.85f
#define MINSNR_F 0.001f
#define MAXSNR_F 1000.0f
#define ONEMINUS_SNR_SMOOTHER_NS_F (1.0f - NS_SNR_SMOOTHER_NS_F)
#define ONEMINUS_MINGAIN_F (1.0f - NS_MINGAIN_F)
#define ONEMINUS_MINGAINU_F (1.0f - NS_MINGAINU_F)
#define ONEMINUS_BETA_F (1.0f - NS_BETA_F)
	class MsrNS
	{
	public:
		MsrNS();
		~MsrNS();
		void process(float* fft_ptr, MsrVAD* pVAD, bool enableNS);
	private:

		float piPriorSNR[FRAME_SIZE];
		float piGain[FRAME_SIZE];
		// Scratch space
		float piPosteriorSNR_NS[FRAME_SIZE];
		float piMLPriorSNR_NS[FRAME_SIZE];
		float piPriorSNR_NS[FRAME_SIZE];
	};
}

#endif /* MSRNS_H_ */
