#include "Pipeline.h"

namespace Beam{
	Pipeline* Pipeline::p_instance = NULL;

	Pipeline::Pipeline() : m_noise_floor(20.0, 0.04, 30000.0, 0.0){
		// initialize noise suppressors.
		for (int channel = 0; channel < MAX_MICROPHONES; ++channel){
			m_ssl_noise_suppressor[channel].init(SAMPLE_RATE, FRAME_SIZE, 1.f, 10.f);
			m_pre_noise_suppressor[channel].init(SAMPLE_RATE, FRAME_SIZE, 1.f, 1.f);
			m_out_noise_suppressor[channel].init(SAMPLE_RATE, FRAME_SIZE, 1.f, 1.f);
		}
		DSPFilter::band_pass_mclt(m_ssl_band_pass_filter, 500.0 / SAMPLE_RATE, 1000.0 / SAMPLE_RATE, 2000.0 / SAMPLE_RATE, 3500.0 / SAMPLE_RATE);
		m_ssl.init(SAMPLE_RATE, FRAME_SIZE);
		// initialize m_frequency_filter.
		for (int sub = 0; sub < MAX_GAIN_SUBBANDS; ++sub){
			DSPFilter::band_pass_mclt(m_frequency_filter[sub], KinectConfig::frequency_bands[sub][0] / SAMPLE_RATE, KinectConfig::frequency_bands[sub][1] / SAMPLE_RATE, KinectConfig::frequency_bands[sub][1] / SAMPLE_RATE, KinectConfig::frequency_bands[sub][2] / SAMPLE_RATE);
		}
		// initialize m_working_frequency.
		m_working_frequency.assign(FRAME_SIZE, std::complex<float>(0.f, 0.f));
		// initialize persistent and dynamic gains
		for (int channel = 0; channel < MAX_MICROPHONES; ++channel){
			m_dynamic_gains[channel].assign(FRAME_SIZE, std::complex<float>(1.f, 0.f));
		}
		for (int channel = 0; channel < MAX_MICROPHONES; ++channel){
			for (int sub = 0; sub < MAX_GAIN_SUBBANDS; ++sub){
				m_persistent_gains[channel][sub] = std::complex<float>(1.f, 0.f);
			}
		}
		m_refresh_ini = 0;
		m_confidence = 0.f;

		// initialize the first and last bin
		m_first_bin = KinectConfig::kinect_descriptor.freq_low / (float)SAMPLE_RATE * (float)FRAME_SIZE * 2.f;
		m_last_bin = KinectConfig::kinect_descriptor.freq_high / (float)SAMPLE_RATE * (float)FRAME_SIZE * 2.f;
	}

	Pipeline* Pipeline::instance(){
		if (p_instance == NULL){
			p_instance = new Pipeline;
		}
		return p_instance;
	}

	void Pipeline::load_profile(){
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

	void Pipeline::preprocess(std::vector<std::complex<float> >* input){
		for (int channel = 0; channel < MAX_MICROPHONES; ++channel){
			for (int bin = 0; bin < FRAME_SIZE; ++bin){
				input[channel][bin] *= m_dynamic_gains[channel][bin];
			}
			m_pre_noise_suppressor[channel].phase_compensation(input[channel]);
		}
	}

	bool Pipeline::source_localize(std::vector<std::complex<float> >* input, double time, float* p_angle){
		//  Apply the SSL band pass filter to the input channels
		//  and have a separate copy of the input channels 
		//  for SSL purposes only
		for (int channel = 0; channel < MAX_MICROPHONES; ++channel){
			for (size_t bin = 0; bin < m_ssl_band_pass_filter.size(); ++bin){
				input[channel][bin] *= m_ssl_band_pass_filter[bin];
			}
		}
		//  Noise suppression
		//  We do heavy noise suppression as we don't care about the musical noises
		//  but we do cary to suppress stationaty noises
		double energy = 0.0;
		for (int channel = 0; channel < MAX_MICROPHONES; ++channel){
			m_ssl_noise_suppressor[channel].noise_compensation(input[channel]);
			energy += (double)Utils::computeRMS(input[channel]);
		}
		energy /= MAX_MICROPHONES;
		double floor = m_noise_floor.nextLevel(time, energy);
		if (energy > 5.290792 * floor){	// TODO: modify this threshold
			// sound signal
			float angle;
			float weight;
			float confidence;
			float std_dev;
			int valid;
			int num;
			m_ssl.process(input, &angle, &weight);
			std::cout << "weight: " << weight << std::endl;
			if (weight > 5e-6f){ // TODO: modify this threshold
				m_ssl.process_next_sample(time, angle, weight);
				m_ssl.get_average(time, p_angle, &confidence, &std_dev, &num, &valid);
				return true;
			}
		}
		return false;
	}

	void Pipeline::beamformer(std::vector<std::complex<float> >* input, std::vector<std::complex<float> >& output, double time){
		//  we have sound source detected - single beam mode
		if (m_confidence > SSL_BEAMCHANGE_CONFIDENCE_THRESHOLD){
			// find the best beam
			float min_dist = FLT_MAX;
			int beam = 0;
			for (int index = 0; index < MAX_BEAMS; ++index){
				float dist = KinectConfig::kinect_beams[index].fi - m_source_position;
				while (dist >(float)PI) dist -= TWO_PI;
				while (dist <= (float)-PI) dist += TWO_PI;
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
				float dist = fabs(KinectConfig::kinect_beams[m_beam].fi - m_source_position);
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
		for (int bin = 0; bin < FRAME_SIZE; ++bin){
			output[bin] = std::complex<float>(0.f, 0.f);
		}
		// TODO: check the hard coding..
		int beg_bin = 6;
		int end_bin = 225;
		for (int bin = beg_bin; bin < end_bin; ++bin){
			std::complex<float> wo0, wo1, wo2, wo3;
			std::complex<float> m0, m1, m2, m3;
			std::complex<float> w0, w1, w2, w3;
			std::complex<float> scale;
			int weight_index = m_beam * FRAME_SIZE * MAX_MICROPHONES + bin * MAX_MICROPHONES;

			m0 = input[0][bin] * scale;
			m1 = input[1][bin] * scale;
			w0 = (KinectConfig::kinect_weights.dd + weight_index)[0];
			w1 = (KinectConfig::kinect_weights.dd + weight_index)[1];
			wo0 = m_pcm_weights[m_beam][0][bin];
			wo1 = m_pcm_weights[m_beam][1][bin];
			m2 = input[2][bin] * scale;
			m3 = input[3][bin] * scale;
			w2 = (KinectConfig::kinect_weights.dd + weight_index)[2];
			w3 = (KinectConfig::kinect_weights.dd + weight_index)[3];
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

	float Pipeline::smart_calibration(float sound_source, std::vector<std::complex<float> >* input){
		float channel_rms[MAX_MICROPHONES] = { 0.f };
		float est_channel_rms[MAX_MICROPHONES] = { 0.f };
		float est_gains[MAX_MICROPHONES] = { 0.f };
		float average_gain = 0.f;
		float sigma = -1.f;
		float theta_rad = 0.f;
		RCoords mic;
		//  Project microphones to the line pointing to the sound source
		//  Here we assume flat wave propagation from the sound source
		for (int channel = 0; channel < MAX_MICROPHONES; ++channel){
			Utils::c2r(mic, KinectConfig::kinect_descriptor.mic[channel].x, KinectConfig::kinect_descriptor.mic[channel].y, KinectConfig::kinect_descriptor.mic[channel].z);
			m_coordinates[channel] = mic.rho * cosf(sound_source - mic.fi) * cosf(mic.theta);
		}
		for (int sub = 0; sub < MAX_GAIN_SUBBANDS; ++sub){
			for (int channel = 0; channel < MAX_MICROPHONES; ++channel){
				for (int bin = 0; bin < FRAME_SIZE; ++bin){
					m_working_frequency[bin] = input[channel][bin] * m_frequency_filter[sub][bin];
				}
				channel_rms[channel] = Utils::computeRMS(m_working_frequency);
			}
			sigma = Utils::approx(m_coordinates, channel_rms, 1, m_coeff, MAX_MICROPHONES);
			float average_rms = 0.f;
			for (int channel = 0; channel < MAX_MICROPHONES; ++channel){
				average_rms += channel_rms[channel];
			}
			average_rms /= MAX_MICROPHONES;
			if (average_rms > FLT_MIN){
				sigma /= average_rms;
			}
			if (m_coeff[1] < 0.f){
				m_coeff[1] = 0.f;
				m_coeff[0] = average_rms;
				sigma = -1.f;
			}
			for (int channel = 0; channel < MAX_MICROPHONES; ++channel){
				est_channel_rms[channel] = m_coeff[1] * m_coordinates[channel] + m_coeff[0];
			}
			for (int channel = 0; channel < MAX_MICROPHONES; ++channel){
				est_gains[channel] = 1.f;
				if (channel_rms[channel] > FLT_MIN){
					est_gains[channel] = est_channel_rms[channel] / channel_rms[channel];
				}
			}
			for (int channel = 0; channel < MAX_MICROPHONES; ++channel){
				average_gain += est_gains[channel];
			}
			average_gain /= MAX_MICROPHONES;
			for (int channel = 0; channel < MAX_MICROPHONES; ++channel){
				est_gains[channel] /= average_gain;
			}
			float weight = 0.001f;
			for (int channel = 0; channel < MAX_MICROPHONES; ++channel){
				if (average_rms > FLT_MIN){
					weight = 0.001f * channel_rms[channel] / average_rms;
				}
				else{
					weight = 0.f;
				}
				float re = m_persistent_gains[channel][sub].real();
				float im = m_persistent_gains[channel][sub].imag();
				m_persistent_gains[channel][sub] = std::complex<float>(re + weight * (est_gains[channel] - re), im);
			}
		}
		++m_refresh_ini;
		if (m_refresh_ini > 200){
			expand_gain();
			m_refresh_ini = 0;
		}
		return sigma;
	}

	void Pipeline::expand_gain(){
		std::complex<float> zero(0.f, 0.f);
		// use MCLT
		float freq_step = (float)SAMPLE_RATE / FRAME_SIZE / 2.f;
		float freq_beg = freq_step / 2.f;
		for (int index = 0; index < FRAME_SIZE; ++index){
			int interp_high = 1;
			int interp_low = 0;
			float freq = freq_beg + index * freq_step;
			if (freq > KinectConfig::frequency_bands[MAX_GAIN_SUBBANDS - 1][0]){
				freq = KinectConfig::frequency_bands[MAX_GAIN_SUBBANDS - 1][0];
			}
			while (freq >= KinectConfig::frequency_bands[interp_high][0] && interp_high < MAX_GAIN_SUBBANDS - 1){
				++interp_high;
			}
			interp_low = interp_high - 1;
			float t = (freq - KinectConfig::frequency_bands[interp_low][0]) / (KinectConfig::frequency_bands[interp_high][0] - KinectConfig::frequency_bands[interp_low][0]);
			// special case 1.  Frequency is less than frequency_bands[0][0] - linear interpolate between 0 & frequency_bands[0][0]
			if (freq < KinectConfig::frequency_bands[0][0]){
				t = freq / KinectConfig::frequency_bands[0][0];
				for (int channel = 0; channel < MAX_MICROPHONES; ++channel){
					m_dynamic_gains[channel][index] = Utils::interpolate(zero, m_persistent_gains[channel][0], t);
				}
			}
			// special case 2, no need to interpolate
			else if (t == 0.f || t >= 1.f){
				int freq_index = t > 0.f ? interp_high : interp_low;
				for (int channel = 0; channel < MAX_MICROPHONES; ++channel){
					m_dynamic_gains[channel][index] = m_persistent_gains[channel][freq_index];
				}
			}
			// standard case  | here we need to interpolate the values
			else{
				for (int channel = 0; channel < MAX_MICROPHONES; ++channel){
					m_dynamic_gains[channel][index] = Utils::interpolate(m_persistent_gains[channel][interp_low], m_persistent_gains[channel][interp_high], t);
				}
			}
		}
	}

	void Pipeline::ansi_bf_msr_process_quad_loop_fast(std::complex<float>* wo0, std::complex<float>* wo1, std::complex<float>* wo2, std::complex<float>* wo3, std::complex<float>& m0, std::complex<float>& m1, std::complex<float>& m2, std::complex<float>& m3, std::complex<float>& w0, std::complex<float>& w1, std::complex<float>& w2, std::complex<float>& w3, float nu, float mu){
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
