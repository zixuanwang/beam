#ifndef DEREVERB_H_
#define DEREVERB_H_

#include "GlobalConfig.h"
#include "Utils.h"
#include <list>

#define TAIL_FRAME_SIZE 3

namespace Beam{
	class DeReverb {
	public:
		DeReverb();
		~DeReverb();
		// cepstral mean subtraction.
		void normalize_cepstral(std::vector<std::complex<float> >& input, bool voice_found);
		// reverbration suppression.
		void suppress(std::vector<std::complex<float> >& input);
		// compute tau for one frequency bin.
		float compute_tau(float init_energy, const std::vector<float>& tail_energy);
		void update_tau(const std::vector<float>& tau);
	private:
		std::vector<std::complex<float> > m_cepstral_mean;
		int m_voice_frame_count;
		bool m_voice_found;
		bool m_tail_found;
		int m_tail_count;
		std::vector<float> m_init_energy;
		std::vector<float> m_energy[FRAME_SIZE];
		std::vector<float> m_tau;
		std::vector<std::list<float> > m_energy_list;
	};
}

#endif /* DEREVERB_H_ */
