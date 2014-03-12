#ifndef MATH_H_
#define MATH_H_

#include <vector>
#include <unordered_map>

class Math {
public:
	Math();
	~Math();

	/// compute l1 norm of the input vector.
	static float l1_norm(const std::vector<float>& array);
	/// l1 normalization.
	static void l1_normalize(std::vector<float>& array);
	/// compute l2 norm of the input vector.
	static float l2_norm(const std::vector<float>& array);
	/// l2 normalization.
	static void l2_normalize(std::vector<float>& array);
	/// normalize angle to [0, 2pi).
	static float normalize_angle(float angle);
	/// square root for each element
	static void square_root(std::vector<float>& array);
	// interpolate the qualtratic function.
	static float interpolateMax(const std::vector<float>& x, const std::vector<float>& y);
	const static float PI;
	const static float HALF_PI;
	const static float TWO_PI;
	const static float TO_RAD;
};

#endif /* MATH_H_ */
