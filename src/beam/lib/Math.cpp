/*
 * Math.cpp
 *
 *  Created on: May 30, 2013
 *      Author: daniewang
 */

#include "Math.h"

const float Math::PI = 3.1415926f;
const float Math::HALF_PI = 1.5707963f;
const float Math::TWO_PI = 6.2831853f;
const float Math::TO_RAD = 0.017453292f;

Math::Math() {

}

Math::~Math() {

}

float Math::l1_norm(const std::vector<float>& array) {
	float sum = 0.f;
	for (auto v : array){
		sum += abs(v);
	}
	return sum;
}

void Math::l1_normalize(std::vector<float>& array) {
	float norm = l1_norm(array);
	if (norm == 0.f)
		return;
	for (auto& v : array){
		v /= norm;
	}
}

float Math::l2_norm(const std::vector<float>& array) {
	float sum = 0.f;
	for (auto v : array){
		sum += v * v;
	}
	return sqrtf(sum);
}

void Math::l2_normalize(std::vector<float>& array) {
	float norm = l2_norm(array);
	if (norm == 0.f)
		return;
	for (auto& v : array){
		v /= norm;
	}
}

void Math::square_root(std::vector<float>& array) {
	for (auto& v : array){
		v = sqrtf(v);
	}
}

float Math::normalize_angle(float angle){
	while (angle < 0.0f) angle += TWO_PI;
	while (angle >= TWO_PI) angle -= TWO_PI;
	return angle;
}

float Math::interpolateMax(const std::vector<float>& x, const std::vector<float>& y){
	float dY20;
	float dY10;
	float dA;
	float dB;
	float dXmax;
	//  Y = A.X2 + B.X + C
	dY20 = (y[2] - y[0]) / (x[2] - x[0]);
	dY10 = (y[1] - y[0]) / (x[1] - x[0]);
	dA = (dY20 - dY10) / (x[2] - x[1]);
	dB = dY20 - dA * (x[2] + x[0]);
	if (dA >= 0.f){
		//  not a maximum??
		return x[1];
	}
	//  the maximum is the zero of the first derivative
	dXmax = -dB / (2.f * dA);
	//  check the result
	if (dXmax > x[2]){
		dXmax = x[1];
	}
	if (dXmax < x[0]){
		dXmax = x[1];
	}
	return dXmax;
}