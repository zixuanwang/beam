#include "Beamformer.h"

namespace Beam{
	Beamformer::Beamformer() : m_beam(0){
		// initialize the first and last bin
		m_first_bin = (int)(KinectConfig::kinect_descriptor.freq_low / (float)SAMPLE_RATE * (float)FRAME_SIZE * 2.f);
		m_last_bin = (int)(KinectConfig::kinect_descriptor.freq_high / (float)SAMPLE_RATE * (float)FRAME_SIZE * 2.f);
	}

	Beamformer::~Beamformer(){
	
	}

	void Beamformer::init(){
		for (int beam = 0; beam < MAX_BEAMS; ++beam){
			for (int channel = 0; channel < MAX_MICROPHONES; ++channel){
				m_pcm_weights[beam][channel].assign(FRAME_SIZE, std::complex<float>(0.f, 0.f));
			}
		}
		// initialize kinect weights
		std::complex<float> zero(0.f, 0.f);
		float freq_step = (float)SAMPLE_RATE / FRAME_SIZE / 2.f;
		float freq_beg = freq_step / 2.f;
		for (int beam = 0; beam < MAX_BEAMS; ++beam){
			int interp_low = 0;
			int interp_high = 1;
			while (KinectConfig::kinect_descriptor.freq_low >= KinectConfig::kinect_weights.frequencies[interp_low] && interp_high < KinectConfig::kinect_weights.num_frequency_bins - 1){
				++interp_high;
				++interp_low;
			}
			for (int bin = 0; bin < FRAME_SIZE; ++bin){
				float freq = freq_beg + bin * freq_step;
				if (freq > KinectConfig::kinect_descriptor.freq_high){
					freq = KinectConfig::kinect_descriptor.freq_high;
				}
				while (freq >= KinectConfig::kinect_weights.frequencies[interp_high] && interp_high < KinectConfig::kinect_weights.num_frequency_bins - 1){
					++interp_high;
				}
				interp_low = interp_high - 1;
				float t = (freq - KinectConfig::kinect_weights.frequencies[interp_low]) / (KinectConfig::kinect_weights.frequencies[interp_high] - KinectConfig::kinect_weights.frequencies[interp_low]);
				// special case 1.  Frequency is less than dFreq_Lo --- set value to 0
				if (freq <= KinectConfig::kinect_descriptor.freq_low){
					for (int channel = 0; channel < MAX_MICROPHONES; ++channel){
						m_pcm_weights[beam][channel][bin] = zero;
					}
				}
				// Special Case 2 - interpolate between 0 and freqency
				else if (t < 0.f){
					int freq_index = KinectConfig::kinect_weights.frequencies[interp_low] > KinectConfig::kinect_descriptor.freq_low ? interp_low : interp_high;
					t = (freq - KinectConfig::kinect_descriptor.freq_low) / (KinectConfig::kinect_weights.frequencies[freq_index] - KinectConfig::kinect_descriptor.freq_low);
					for (int channel = 0; channel < MAX_MICROPHONES; ++channel){
						int weight_index = beam * FRAME_SIZE * MAX_MICROPHONES;
						weight_index += freq_index * KinectConfig::kinect_weights.num_channels;
						weight_index += channel;
						m_pcm_weights[beam][channel][bin] = Utils::interpolate(zero, KinectConfig::kinect_weights.weights[weight_index + KinectConfig::kinect_weights.num_channels], t);
					}
				}
				// special case 3, no need to interpolate
				else if (t == 0.f || t >= 1.f){
					int freq_index = t > 0.f ? interp_high : interp_low;
					for (int channel = 0; channel < MAX_MICROPHONES; ++channel){
						int weight_index = beam * KinectConfig::kinect_weights.num_frequency_bins * KinectConfig::kinect_weights.num_channels;
						weight_index += freq_index * KinectConfig::kinect_weights.num_channels;
						weight_index += channel;
						m_pcm_weights[beam][channel][bin] = KinectConfig::kinect_weights.weights[weight_index];
					}
				}
				// standard case  | here we need to interpolate the values
				else{
					for (int channel = 0; channel < MAX_MICROPHONES; ++channel){
						int weight_index = beam * KinectConfig::kinect_weights.num_frequency_bins * KinectConfig::kinect_weights.num_channels;
						weight_index += interp_low * KinectConfig::kinect_weights.num_channels;
						weight_index += channel;
						m_pcm_weights[beam][channel][bin] = Utils::interpolate(KinectConfig::kinect_weights.weights[weight_index], KinectConfig::kinect_weights.weights[weight_index + KinectConfig::kinect_weights.num_channels], t);
					}
				}
			}
		}
	}

	void Beamformer::compute(std::vector<std::complex<float> >* input, std::vector<std::complex<float> >& output, float angle, double time){
		//  we have sound source detected - single beam mode
		// find the best beam
		float min_dist = FLT_MAX;
		int beam = 0;
		for (int index = 0; index < MAX_BEAMS; ++index){
			float dist = KinectConfig::kinect_beams[index].fi - angle;
			while (dist >(float)PI) dist -= (float)TWO_PI;
			while (dist <= (float)-PI) dist += (float)TWO_PI;
			dist = fabs(dist);
			if (dist < min_dist){
				min_dist = dist;
				beam = index;
			}
		}
		//  big difference - switch the beam instantly
		if (abs(m_beam - beam) > 1){
			m_beam = beam;
		}
		else{
			//  neighbor beams - switch only if two thirds of the way
			float dist = fabs(KinectConfig::kinect_beams[m_beam].fi - angle);
			if (beam == 0){
				if (dist > 0.66f * (KinectConfig::kinect_beams[beam + 1].fi - KinectConfig::kinect_beams[beam].fi)){
					m_beam = beam;
				}
			}
			else if (beam == MAX_BEAMS - 1){
				if (dist > 0.66f * (KinectConfig::kinect_beams[beam].fi - KinectConfig::kinect_beams[beam - 1].fi)){
					m_beam = beam;
				}
			}
			else{
				if (dist > 0.66f * (KinectConfig::kinect_beams[beam + 1].fi - KinectConfig::kinect_beams[beam].fi) || dist > 0.66f * (KinectConfig::kinect_beams[beam].fi - KinectConfig::kinect_beams[beam - 1].fi)){
					m_beam = beam;
				}
			}
		}
		for (int bin = 0; bin < m_first_bin; ++bin){
			output[bin] = std::complex<float>(0.f, 0.f);
		}
		for (int bin = m_first_bin; bin < m_last_bin; ++bin){
			output[bin] = std::complex<float>(0.f, 0.f);
			for (int channel = 0; channel < MAX_MICROPHONES; ++channel){
				output[bin] += m_pcm_weights[m_beam][channel][bin] * input[channel][bin];
			}
		}
		for (int bin = m_last_bin; bin < FRAME_SIZE; ++bin){
			output[bin] = std::complex<float>(0.f, 0.f);
		}
		int beg_bin = 6;
		int end_bin = 225;
		for (int bin = beg_bin; bin < end_bin; ++bin){
			std::complex<float> wo0, wo1, wo2, wo3;
			std::complex<float> m0, m1, m2, m3;
			std::complex<float> w0, w1, w2, w3;
			std::complex<float> scale(F2RAISED23_INV, 0.f);
			int weight_index = m_beam * FRAME_SIZE * MAX_MICROPHONES + bin * MAX_MICROPHONES;

			m0 = input[0][bin] * scale;
			m1 = input[1][bin] * scale;
			m2 = input[2][bin] * scale;
			m3 = input[3][bin] * scale;
			w0 = (KinectConfig::kinect_weights.dd + weight_index)[0];
			w1 = (KinectConfig::kinect_weights.dd + weight_index)[1];
			w2 = (KinectConfig::kinect_weights.dd + weight_index)[2];
			w3 = (KinectConfig::kinect_weights.dd + weight_index)[3];
			wo0 = m_pcm_weights[m_beam][0][bin];
			wo1 = m_pcm_weights[m_beam][1][bin];
			wo2 = m_pcm_weights[m_beam][2][bin];
			wo3 = m_pcm_weights[m_beam][3][bin];
			ansi_bf_msr_process_quad_loop_fast(&wo0, &wo1, &wo2, &wo3, m0, m1, m2, m3, w0, w1, w2, w3, 0.0000010f, 0.0080000f);
			// Update working weights
			m_pcm_weights[m_beam][0][bin] = wo0;
			m_pcm_weights[m_beam][1][bin] = wo1;
			m_pcm_weights[m_beam][2][bin] = wo2;
			m_pcm_weights[m_beam][3][bin] = wo3;
		}
	}

	void Beamformer::ansi_bf_msr_process_quad_loop_fast(std::complex<float>* wo0, std::complex<float>* wo1, std::complex<float>* wo2, std::complex<float>* wo3, std::complex<float>& m0, std::complex<float>& m1, std::complex<float>& m2, std::complex<float>& m3, std::complex<float>& w0, std::complex<float>& w1, std::complex<float>& w2, std::complex<float>& w3, float nu, float mu){
		/* Basic idea for this code:                                         */
		/* 1) Pull out the current values for each microphone for this bin   */
		/* 2) Find the complex outer product (x * transpose(conj(x)))        */
		/* 3) Take NumMics by NumMics matrix and add mu to the diagonal      */
		/* 4) Invert the matrix                                              */
		/* 5) Multiply the orignal weights (complex transpose) by the matrix */
		/* 6) Divide by (Orig complex transpose)*(inverted Matric)*(Orig)    */
		/* 7) Finally weight the stored values back into the frequency bin   */
		/* Create the denominator first. */
		std::complex<float> wn0, wn1, wn2, wn3;
		float den = 0.f;
		float tmp = 0.f;
		/* To understand the notation imagine the microphones are:       */
		/* M0 = (A,B); M1 = (C,D); M2 = (E,F); M3 = (G,H)                */
		/* Imagine the weights are:                                      */
		/* W0 = (J,K); W1 = (L,M); W2 = (O,P); W3 = (Q,R)                */

		/* This algorithm was found by writing an algebraic maniuplation */
		/* program from scratch and then analyzing the output of the     */
		/* program.  From this common terms are thrown out and the det   */
		/* of the matrix is thrown away since the same det is on the top */
		/* and the bottom.  Additionally all terms are biased by Nu*Nu   */
		/* and so that factor is thrown out as well.                     */

		/* First we will take care of the power component of the den and */
		/* the weights times the diagonal to set-up the output weights.  */
		float out0 = nu;
		float out1 = nu;
		float out2 = nu;
		float out3 = nu;

		tmp = std::norm(m0);
		out1 += tmp;
		out2 += tmp;
		out3 += tmp;
		tmp = std::norm(m1);
		out0 += tmp;
		out2 += tmp;
		out3 += tmp;
		tmp = std::norm(m2);
		out0 += tmp;
		out1 += tmp;
		out3 += tmp;
		tmp = std::norm(m3);
		out0 += tmp;
		out1 += tmp;
		out2 += tmp;
		den += out0 * std::norm(w0);
		den += out1 * std::norm(w1);
		den += out2 * std::norm(w2);
		den += out3 * std::norm(w3);
		wn0 = std::complex<float>(out0 * w0.real(), -1.f * out0 * w0.imag());
		wn1 = std::complex<float>(out1 * w1.real(), -1.f * out1 * w1.imag());
		wn2 = std::complex<float>(out2 * w2.real(), -1.f * out2 * w2.imag());
		wn3 = std::complex<float>(out3 * w3.real(), -1.f * out3 * w3.imag());
		/* Need to put in the mixed terms now. */
		// AC & BD
		tmp = -1.f * m0.real() * m1.real() - m0.imag() * m1.imag();
		wn0 += std::complex<float>(tmp * w1.real(), -tmp * w1.imag());
		wn1 += std::complex<float>(tmp * w0.real(), -tmp * w0.imag());
		den += 2.f * tmp * (w0.real() * w1.real() + w0.imag() * w1.imag());
		// BC & AD
		tmp = m0.imag() * m1.real() - m0.real() * m1.imag();
		wn0 += std::complex<float>(tmp * w1.imag(), tmp * w1.real());
		wn1 -= std::complex<float>(tmp * w0.imag(), tmp * w0.real());
		den += 2.f * tmp * (w0.real() * w1.imag() - w0.imag() * w1.real());
		// BF & AE
		tmp = -1.f * m0.real() * m2.real() - m0.imag() * m2.imag();
		wn0 += std::complex<float>(tmp * w2.real(), -tmp * w2.imag());
		wn2 += std::complex<float>(tmp * w0.real(), -tmp * w0.imag());
		den += 2.f * tmp * (w0.real() * w2.real() + w0.imag() * w2.imag());
		// BE & AF
		tmp = m0.imag() * m2.real() - m0.real() * m2.imag();
		wn0 += std::complex<float>(tmp * w2.imag(), tmp * w2.real());
		wn2 -= std::complex<float>(tmp * w0.imag(), tmp * w0.real());
		den += 2.f * tmp * (w0.real() * w2.imag() - w0.imag() * w2.real());
		// BH & AG
		tmp = -1.f * m0.real() * m3.real() - m0.imag() * m3.imag();
		wn0 += std::complex<float>(tmp * w3.real(), -tmp * w3.imag());
		wn3 += std::complex<float>(tmp * w0.real(), -tmp * w0.imag());
		den += 2.f * tmp * (w0.real() * w3.real() + w0.imag() * w3.imag());
		// BG & AH
		tmp = m0.imag() * m3.real() - m0.real() * m3.imag();
		wn0 += std::complex<float>(tmp * w3.imag(), tmp * w3.real());
		wn3 -= std::complex<float>(tmp * w0.imag(), tmp * w0.real());
		den += 2.f * tmp * (w0.real() * w3.imag() - w0.imag() * w3.real());
		// DF & CE
		tmp = -1.f * m1.real() * m2.real() - m1.imag() * m2.imag();
		wn1 += std::complex<float>(tmp * w2.real(), -tmp * w2.imag());
		wn2 += std::complex<float>(tmp * w1.real(), -tmp * w1.imag());
		den += 2.f * tmp * (w1.real() * w2.real() + w1.imag() * w2.imag());
		// DE & CF
		tmp = m1.imag() * m2.real() - m1.real() * m2.imag();
		wn1 += std::complex<float>(tmp * w2.imag(), tmp * w2.real());
		wn2 -= std::complex<float>(tmp * w1.imag(), tmp * w1.real());
		den += 2.f * tmp * (w1.real() * w2.imag() - w1.imag() * w2.real());
		// CG & DH
		tmp = -1.f * m1.real() * m3.real() - m1.imag() * m3.imag();
		wn1 += std::complex<float>(tmp * w3.real(), -tmp * w3.imag());
		wn3 += std::complex<float>(tmp * w1.real(), -tmp * w1.imag());
		den += 2.f * tmp * (w1.real() * w3.real() + w1.imag() * w3.imag());
		// CH & DG
		tmp = m1.imag() * m3.real() - m1.real() * m3.imag();
		wn1 += std::complex<float>(tmp * w3.imag(), tmp * w3.real());
		wn3 -= std::complex<float>(tmp * w1.imag(), tmp * w1.real());
		den += 2.f * tmp * (w1.real() * w3.imag() - w1.imag() * w3.real());
		// EG & FH
		tmp = -1.f * m2.real() * m3.real() - m2.imag() * m3.imag();
		wn2 += std::complex<float>(tmp * w3.real(), -tmp * w3.imag());
		wn3 += std::complex<float>(tmp * w2.real(), -tmp * w2.imag());
		den += 2.f * tmp * (w2.real() * w3.real() + w2.imag() * w3.real());
		// EH & FG
		tmp = m2.imag() * m3.real() - m2.real() * m3.imag();
		wn2 += std::complex<float>(tmp * w3.imag(), tmp * w3.real());
		wn3 -= std::complex<float>(tmp * w2.imag(), tmp * w2.real());
		den += 2.f * tmp * (w2.real() * w3.imag() - w2.imag() * w3.real());
		/* Now divide by the denominator. */
		wn0 /= den;
		wn1 /= den;
		wn2 /= den;
		wn3 /= den;
		/* So now we have the new weights so now adjust the weights for the beamformer. */
		*wo0 = *wo0 * (1.f - mu) + mu * wn0;
		*wo1 = *wo1 * (1.f - mu) + mu * wn1;
		*wo2 = *wo2 * (1.f - mu) + mu * wn2;
		*wo3 = *wo3 * (1.f - mu) + mu * wn3;
	}
}