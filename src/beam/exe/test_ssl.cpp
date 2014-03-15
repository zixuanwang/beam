#include <iostream>
#include "beam\lib\MCLT.h"
#include "beam\lib\SoundSourceLocalizer.h"


int main(int argc, char* argv[]){
	int M = 4;
	MCLT mclt(M);
	double input[8];
	input[0] = 1.0;
	input[1] = 2.0;
	input[2] = 3.0;
	input[3] = 4.0;
	input[4] = 6.0;
	input[5] = 8.0;
	input[6] = -1.0;
	input[7] = 23.0;
	double re[4];
	double im[4];
	mclt.compute(input, re, im);

	return 0;
}