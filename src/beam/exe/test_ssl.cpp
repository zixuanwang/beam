#include <iostream>
#include "beam\lib\kiss_fft.h"
#include "beam\lib\SoundSourceLocalizer.h"


int main(int argc, char* argv[]){
	kiss_fft_cfg cfg = kiss_fft_alloc(5, 0, 0, 0);
	kiss_fft_cpx input[5] = { kiss_fft_cpx{ 1.f, 0.f }, kiss_fft_cpx{ 2.f, 0.f }, kiss_fft_cpx{ 3.f, 0.f }, kiss_fft_cpx{ 4.f, 0.f }, kiss_fft_cpx{ 5.f, 0.f } };
	kiss_fft_cpx output[5];
	kiss_fft(cfg, input, output);
	free(cfg);
	return 0;
}