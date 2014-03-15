#include <iostream>
#include "beam\lib\MCLT.h"
#include "beam\lib\SoundSourceLocalizer.h"


int main(int argc, char* argv[]){
	int M = 2;
	MCLT mclt(M);
	double input[4];
	input[0] = 1.0;
	input[1] = 2.0;
	input[2] = 3.0;
	input[3] = 4.0;
	std::complex<double> output[2];
	mclt.compute(input, output);
	//mclt.compute(input, output);

	return 0;
}