#include "GSCBeamformer.h"

namespace Beam{
	GSCBeamformer::GSCBeamformer(){
		for (int i = 0; i < MAX_MICROPHONES; ++i){
			std::fill(m_input_prev[i], m_input_prev[i] + FRAME_SIZE, 0.f);
			std::fill(m_y_prev[i], m_y_prev[i] + FRAME_SIZE, 0.f);
			std::fill(m_bm[i], m_bm[i] + BM_N, 1.f / BM_N);
			std::fill(m_mc[i], m_mc[i] + MC_L, 1.f / MC_L);
		}
		std::fill(m_d_prev, m_d_prev + FRAME_SIZE, 0.f);
	}

	GSCBeamformer::~GSCBeamformer(){
	
	}

	void GSCBeamformer::compute(float output[FRAME_SIZE], float input[][FRAME_SIZE], float angle, bool voice, float ref[FRAME_SIZE]){
		// skip delay sum beamformer now.
		// process in the time domain.
		int bm_p = 4;
		int bm_n = 8;
		float d[FRAME_SIZE] = { 0.f };
		float y[MAX_MICROPHONES][FRAME_SIZE] = { 0.f };
		std::fill(output, output + FRAME_SIZE, 0.f);
		if (ref != NULL){
			std::copy(ref, ref + FRAME_SIZE, d);
		}
		else{
			for (int bin = 0; bin < FRAME_SIZE; ++bin){
				for (int channel = 0; channel < MAX_MICROPHONES; ++channel){
					d[bin] += input[channel][bin];
				}
				d[bin] /= MAX_MICROPHONES;
			}
		}
		// gsc bm

		for (int k = 0; k < FRAME_SIZE; ++k){
			if (voice){
				for (int channel = 0; channel < MAX_MICROPHONES; ++channel){
					float x = 0.f;
					if (k - BM_P < 0){
						x = m_input_prev[channel][FRAME_SIZE + k - BM_P];
					}
					else{
						x = input[channel][k - BM_P];
					}
					float sum = 0.f;
					float norm = 0.f;
					for (int j = 0; j < BM_N; ++j){
						if (k - j < 0){
							float v = m_d_prev[FRAME_SIZE + k - j];
							sum += m_bm[channel][j] * v;
							norm += v * v;
						}
						else{
							float v = d[k - j];
							sum += m_bm[channel][j] * v;
							norm += v * v;
						}
					}
					norm = sqrtf(norm);
					float e = x - sum;
					if (norm != 0.f){
						for (int j = 0; j < BM_N; ++j){
							if (k - j < 0){
								m_bm[channel][j] += e / norm * m_d_prev[FRAME_SIZE + k - j];
							}
							else{
								m_bm[channel][j] += e / norm * d[k - j];
							}
						}
					}
				}
			}
			// compute y
			for (int channel = 0; channel < MAX_MICROPHONES; ++channel){
				float x = 0.f;
				if (k - BM_P < 0){
					x = m_input_prev[channel][FRAME_SIZE + k - BM_P];
				}
				else{
					x = input[channel][k - BM_P];
				}
				float sum = 0.f;
				for (int j = 0; j < BM_N; ++j){
					if (k - j < 0){
						float v = m_d_prev[FRAME_SIZE + k - j];
						sum += m_bm[channel][j] * v;
					}
					else{
						float v = d[k - j];
						sum += m_bm[channel][j] * v;
					}
				}
				y[channel][k] = x - sum;
			}
			// gsc mc
			if (!voice){
				float z = 0.f;
				float norm = 0.f;
				for (int channel = 0; channel < MAX_MICROPHONES; ++channel){
					for (int j = 0; j < MC_L; ++j){
						if (k - j < 0){
							z -= m_y_prev[channel][FRAME_SIZE + k - j] * m_mc[channel][j];
							norm += m_y_prev[channel][FRAME_SIZE + k - j] * m_y_prev[channel][FRAME_SIZE + k - j];
						}
						else{
							z -= y[channel][k - j] * m_mc[channel][j];
							norm += y[channel][k - j] * y[channel][k - j];
						}
					}
				}
				if (k - MC_Q < 0){
					z += m_d_prev[FRAME_SIZE + k - MC_Q];
				}
				else{
					z += d[k - MC_Q];
				}
				if (norm != 0.f){
					for (int channel = 0; channel < MAX_MICROPHONES; ++channel){
						for (int j = 0; j < MC_L; ++j){
							if (k - j < 0){
								m_mc[channel][j] += z / norm * m_y_prev[channel][FRAME_SIZE + k - j];
							}
							else{
								m_mc[channel][j] += z / norm * y[channel][k - j];
							}
						}
					}
				}
			}
			float z = 0.f;
			for (int channel = 0; channel < MAX_MICROPHONES; ++channel){
				for (int j = 0; j < MC_L; ++j){
					if (k - j < 0){
						z -= m_y_prev[channel][FRAME_SIZE + k - j] * m_mc[channel][j];
					}
					else{
						z -= y[channel][k - j] * m_mc[channel][j];
					}
				}
			}
			if (k - MC_Q < 0){
				z += m_d_prev[FRAME_SIZE + k - MC_Q];
			}
			else{
				z += d[k - MC_Q];
			}
			output[k] = z;
		}
		// copy input. copy y.
		for (int channel = 0; channel < MAX_MICROPHONES; ++channel){
			std::copy(input[channel], input[channel] + FRAME_SIZE, m_input_prev[channel]);
			std::copy(y[channel], y[channel] + FRAME_SIZE, m_y_prev[channel]);
		}
		// copy d.
		std::copy(d, d + FRAME_SIZE, m_d_prev);
	}
}