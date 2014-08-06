#include "Pipeline.h"

namespace Beam{
	Pipeline* Pipeline::p_instance = NULL;

	Pipeline::Pipeline() : m_noise_floor(20.0, 0.04, 30000.0, 0.0){
		// initialize band pass filter.
		DSPFilter::band_pass_mclt(m_band_pass_filter, 500.f / SAMPLE_RATE, 1000.f / SAMPLE_RATE, 2000.f / SAMPLE_RATE, 3500.f / SAMPLE_RATE);
		// initialize noise suppressors.
		for (int channel = 0; channel < MAX_MICROPHONES; ++channel){
			m_ssl_noise_suppressor[channel].init(SAMPLE_RATE, FRAME_SIZE, 1.f, 10.f);
			m_pre_noise_suppressor[channel].init(SAMPLE_RATE, FRAME_SIZE, 1.f, 1.f);
			m_pre_suppressor[channel].init(SAMPLE_RATE, FRAME_SIZE, 1.f, 10.f);
		}
		m_out_noise_suppressor.init(SAMPLE_RATE, FRAME_SIZE, 1.f, 10.f);
		m_ssl.init(SAMPLE_RATE, FRAME_SIZE);
		// initialize persistent and dynamic gains
		for (int channel = 0; channel < MAX_MICROPHONES; ++channel){
			m_dynamic_gains[channel].assign(FRAME_SIZE, std::complex<float>(1.f, 0.f));
		}
		for (int channel = 0; channel < MAX_MICROPHONES; ++channel){
			for (int sub = 0; sub < MAX_GAIN_SUBBANDS; ++sub){
				m_persistent_gains[channel][sub] = std::complex<float>(1.f, 0.f);
			}
		}
		// initialize m_input_channels
		for (int channel = 0; channel < MAX_MICROPHONES; ++channel){
			m_input_channels[channel].assign(FRAME_SIZE, std::complex<float>(0.f, 0.f));
		}
		// initialize input buffers
		for (int channel = 0; channel < MAX_MICROPHONES; ++channel){
			std::fill(m_input_prev[channel], m_input_prev[channel] + FRAME_SIZE, 0.f);
			std::fill(m_input[channel], m_input[channel] + 2 * FRAME_SIZE, 0.f);
		}
		std::fill(m_output_prev, m_output_prev + FRAME_SIZE, 0.f);
		std::fill(m_output, m_output + 2 * FRAME_SIZE, 0.f);
		for (int channel = 0; channel < MAX_MICROPHONES; ++channel){
			m_frequency_input[channel].assign(FRAME_SIZE, std::complex<float>(0.f, 0.f));
		}
		m_frequency_output.assign(FRAME_SIZE, std::complex<float>(0.f, 0.f));
		// initialize gains.
		expand_gain();
		m_refresh_gain = 0;
		m_confidence = 0.f;
		// initialize m_time.
		m_time = 0.0;
		// initialize m_angle.
		m_angle = 0.f;
		// initialize m_voice_found.
		m_voice_found = false;
		// initialize m_source_found.
		m_source_found = false;
		// initialize m_frame_number.
		m_frame_number = 0;
	}

	Pipeline* Pipeline::instance(){
		if (p_instance == NULL){
			p_instance = new Pipeline;
		}
		return p_instance;
	}

	void Pipeline::phase_compensation(float* fft_ptr, bool analysis){
		int fftSize = FRAME_SIZE * 2;
		float* realPtr = fft_ptr + 1;
		float* imagPtr = fft_ptr + fftSize - 1;

		if ((((m_frame_number & 0x03) == 1) && analysis) ||
			(((m_frame_number & 0x03) == 3) && !analysis))
		{
			// first sample X[1]
			float fltTemp = *realPtr;
			*realPtr = -*imagPtr;
			*imagPtr = fltTemp;
			realPtr++;
			imagPtr--;
			for (int i = 2; i < fftSize / 2; i = i + 2)
			{
				// index i is corresponding to X[i].
				fltTemp = *realPtr;
				*realPtr = *imagPtr;
				*imagPtr = -fltTemp;
				realPtr++;
				imagPtr--;
				fltTemp = *realPtr;
				*realPtr = -*imagPtr;
				*imagPtr = fltTemp;
				realPtr++;
				imagPtr--;
			}
		}
		else if ((((m_frame_number & 0x03) == 1) && !analysis) ||
			(((m_frame_number & 0x03) == 3) && analysis))
		{
			// First sample, X[1]
			float fltTemp = *realPtr;
			*realPtr = *imagPtr;
			*imagPtr = -fltTemp;
			realPtr++;
			imagPtr--;
			for (int i = 2; i < fftSize / 2; i = i + 2)
			{
				// index i is corresponding to X[i].
				fltTemp = *realPtr;
				*realPtr = -*imagPtr;
				*imagPtr = fltTemp;
				realPtr++;
				imagPtr--;
				fltTemp = *realPtr;
				*realPtr = *imagPtr;
				*imagPtr = -fltTemp;
				realPtr++;
				imagPtr--;
			}
		}
		else if ((m_frame_number & 0x03) == 2)
		{  // phase compensation term is -1
			for (int i = 1; i < fftSize / 2; i++)
			{
				*realPtr = -*realPtr;
				*imagPtr = -*imagPtr;
				realPtr++;
				imagPtr--;
			}
			// CTODO: can not do phase compensation for coefficient fftSize/2 - 1
		}
	}

	void Pipeline::convert_input(std::vector<std::complex<float> >& input, float* fft_ptr){
		int two_frame_size = 2 * FRAME_SIZE;
		input[0].real(fft_ptr[0]);
		input[0].imag(0.f);
		for (int i = 1; i < FRAME_SIZE; ++i){
			input[i].real(fft_ptr[i]);
			input[i].imag(fft_ptr[two_frame_size - i]);
		}
	}

	void Pipeline::convert_output(const std::vector<std::complex<float> >& output, float* fft_ptr){
		int two_frame_size = 2 * FRAME_SIZE;
		fft_ptr[0] = output[0].real();
		fft_ptr[FRAME_SIZE] = 0.f;
		for (int i = 1; i < FRAME_SIZE; ++i){
			fft_ptr[i] = output[i].real();
			fft_ptr[two_frame_size - i] = output[i].imag();
		}
	}

	void Pipeline::process(float input[MAX_MICROPHONES][FRAME_SIZE], float output[FRAME_SIZE]){
		for (int channel = 0; channel < MAX_MICROPHONES; ++channel){
			for (int i = 0; i < FRAME_SIZE; ++i) {
				m_input[channel][i] = m_input_prev[channel][i];
			}
			for (int i = FRAME_SIZE; i < 2 * FRAME_SIZE; ++i) {
				m_input[channel][i] = input[channel][i - FRAME_SIZE];
			}
			std::copy(input[channel], input[channel] + FRAME_SIZE, m_input_prev[channel]);
			float input_fft[2 * FRAME_SIZE];
			MCLT::AecCcsFwdMclt(m_input[channel], input_fft, true);
			Beam::Pipeline::instance()->phase_compensation(input_fft, true);
			Beam::Pipeline::instance()->convert_input(m_frequency_input[channel], input_fft);
		}
		Beam::Pipeline::instance()->preprocess(m_frequency_input); // noise suppression and dynamic gain
		float angle;
		Beam::Pipeline::instance()->source_localize(m_frequency_input, &angle); // sound source localization
		Beam::Pipeline::instance()->smart_calibration(m_frequency_input); // calibration
		Beam::Pipeline::instance()->beamforming(m_frequency_input, m_frequency_output); // beamforming
		Beam::Pipeline::instance()->postprocessing(m_frequency_output); // NS
		float output_fft[2 * FRAME_SIZE];
		Beam::Pipeline::instance()->convert_output(m_frequency_output, output_fft);
		Beam::Pipeline::instance()->phase_compensation(output_fft, false);
		Beam::MCLT::AecCcsInvMclt(output_fft, m_output, true);
		for (int i = 0; i < FRAME_SIZE; ++i){
			output[i] = m_output[i] + m_output_prev[i];
		}
		std::copy(m_output + FRAME_SIZE, m_output + 2 * FRAME_SIZE, m_output_prev);
		++m_frame_number;
	}

	void Pipeline::preprocess(std::vector<std::complex<float> >* input){
		for (int channel = 0; channel < MAX_MICROPHONES; ++channel){
			// TODO check dynamic gains here.
			//m_pre_suppressor[channel].noise_compensation(input[channel]); // NS here.
			for (int bin = 0; bin < FRAME_SIZE; ++bin){
				input[channel][bin] *= m_dynamic_gains[channel][bin];
			}
			m_pre_noise_suppressor[channel].phase_compensation(input[channel]);
		}
	}

	void Pipeline::source_localize(std::vector<std::complex<float> >* input, float* p_angle){
		m_time += (double)FRAME_SIZE / (double)SAMPLE_RATE;
		m_source_found = false;
		//  Apply the SSL band pass filter to the input channels
		//  and have a separate copy of the input channels 
		//  for SSL purposes only
		for (int channel = 0; channel < MAX_MICROPHONES; ++channel){
			for (size_t bin = 0; bin < m_band_pass_filter.size(); ++bin){
				m_input_channels[channel][bin] = input[channel][bin] * m_band_pass_filter[bin];
			}
		}
		//  Noise suppression
		//  We do heavy noise suppression as we don't care about the musical noises
		//  but we do cary to suppress stationaty noises
		double energy = 0.0;
		for (int channel = 0; channel < MAX_MICROPHONES; ++channel){
			//m_ssl_noise_suppressor[channel].noise_compensation(m_input_channels[channel]);
			energy += (double)Utils::computeRMS(m_input_channels[channel]);
		}
		energy /= MAX_MICROPHONES;
		double floor = m_noise_floor.nextLevel(m_time, energy);
		if (energy > SSL_RELATIVE_ENERGY_THRESHOLD * floor && energy > SSL_ABSOLUTE_ENERGY_THRESHOLD){
			// sound signal
			m_voice_found = true;
			float angle;
			float weight;
			m_ssl.process(m_input_channels, input, &angle, &weight);
			if (weight > SSL_CONTRAST_THRESHOLD){
				m_ssl.process_next_sample(m_time, angle, weight);
				m_source_found = true;
			}
		}
		else{
			m_voice_found = false;
		}
		float std_dev;
		int valid;
		int num;
		m_ssl.get_average(m_time, p_angle, &m_confidence, &std_dev, &num, &valid);
		m_angle = *p_angle;
	}

	void Pipeline::dereverbration(std::vector<std::complex<float> >* input){
		for (int channel = 0; channel < MAX_MICROPHONES; ++channel){
			m_dereverb[channel].suppress(input[channel]);
		}
		//m_dereverb[0].suppress(input);
	}

	void Pipeline::smart_calibration(std::vector<std::complex<float> >* input){
		if (m_source_found){
			float sigma = m_calibrator.calibrate(m_angle, m_input_channels, m_persistent_gains);
			++m_refresh_gain;
			if (m_refresh_gain > 200){
				expand_gain();
				m_refresh_gain = 0;
			}
		}
	}

	void Pipeline::beamforming(std::vector<std::complex<float> >* input, std::vector<std::complex<float> >& output){
		m_beamformer.compute(input, output, 0.f, m_confidence, m_time);
	}

	void Pipeline::postprocessing(std::vector<std::complex<float> >& input){
		//m_out_noise_suppressor.frequency_shifting(input);
		m_out_noise_suppressor.noise_compensation(input);
	}

	void Pipeline::expand_gain(){
		std::complex<float> zero(0.f, 0.f);
		// use MCLT
		float freq_step = (float)SAMPLE_RATE / FRAME_SIZE / 2.f;
		float freq_beg = 0.f;
		if (USE_MCLT){
			freq_beg = freq_step / 2.f;
		}
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
}
