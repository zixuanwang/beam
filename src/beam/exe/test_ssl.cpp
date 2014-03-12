#include <iostream>
#include "beam\lib\SoundSourceLocalizer.h"

int main(int argc, char* argv[]){
	float v[5] = { 1.f, 4.f, -1.f, 0.f, 3.f };
	auto p = std::max_element(v, v + 5);
	std::cout << p - v << std::endl;
	return 0;
}