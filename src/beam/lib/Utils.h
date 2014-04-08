#ifndef UTILS_H_
#define UTILS_H_

#include <complex>
#include <memory>
#include <vector>
#include "Coords.h"

namespace Beam{
#define PI 3.1415926535897932384626433832795
#define HALF_PI 1.5707963267948966192313216916398
#define TWO_PI 6.283185307179586476925286766559
#define TO_RAD 0.017453292519943295769236907684886

	class Utils {
	public:
		Utils(){
		}
		~Utils(){
		}
		/// compute l1 norm of the input vector.
		template<typename T>
		static T l1_norm(const std::vector<T>& array){
			T sum = 0;
			for (auto v : array){
				sum += abs(v);
			}
			return sum;
		}
		/// l1 normalization.
		template<typename T>
		static void l1_normalize(std::vector<T>& array){
			T norm = l1_norm(array);
			if (norm == 0)
				return;
			for (auto& v : array){
				v /= norm;
			}
		}
		/// compute l2 norm of the input vector.
		template<typename T>
		static T l2_norm(const std::vector<T>& array){
			T sum = 0;
			for (auto v : array){
				sum += v * v;
			}
			return sqrt(sum);
		}
		/// l2 normalization.
		template<typename T>
		static void l2_normalize(std::vector<T>& array){
			T norm = l2_norm(array);
			if (norm == 0)
				return;
			for (auto& v : array){
				v /= norm;
			}
		}
		/// normalize angle to [0, 2pi).
		template<typename T>
		static T normalize_angle(T angle){
			while (angle < (T)(0)) angle += (T)(TWO_PI);
			while (angle >= (T)(TWO_PI)) angle -= (T)(TWO_PI);
			return angle;
		}
		/// normalize a complex number to have unit amplitude.
		template<typename T>
		static void normalize_complex(std::complex<T>& v){
			T n = std::abs(v);
			if (n >(T)0){
				v /= n;
			}
		}
		/// interpolate two complex numbers
		template<typename T>
		static std::complex<T> interpolate(const std::complex<T>& p0, const std::complex<T>& p1, T t){
			T rho0 = std::abs(p0);
			T rho1 = std::abs(p1);
			T rho = rho0 + (t * (rho1 - rho0));
			T fi0 = normalize_angle(std::arg(p0));
			T fi1 = normalize_angle(std::arg(p1));
			T ffi = fi1 - fi0;
			if (ffi > (T)PI) ffi -= (T)TWO_PI;
			if (ffi <= (T)-PI) ffi += (T)TWO_PI;
			ffi = fi0 + t * ffi;
			std::complex<T> result(rho * cos(ffi), rho * sin(ffi));
			return result;
		}

		template<typename T>
		static T computeRMS(const std::vector<std::complex<T> >& array){
			T energy = (T)0;
			for (auto& v : array){
				energy += std::norm(v);
			}
			return sqrt(energy / (T)array.size());
		}
		/// square root for each element
		template<typename T>
		static void square_root(std::vector<T>& array){
			for (auto& v : array){
				v = sqrt(v);
			}
		}
		// interpolate the qualtratic function.
		template<typename T>
		static T interpolate_max(const std::vector<T>& x, const std::vector<T>& y){
			T dY20;
			T dY10;
			T dA;
			T dB;
			T dXmax;
			//  Y = A.X2 + B.X + C
			dY20 = (y[2] - y[0]) / (x[2] - x[0]);
			dY10 = (y[1] - y[0]) / (x[1] - x[0]);
			dA = (dY20 - dY10) / (x[2] - x[1]);
			dB = dY20 - dA * (x[2] + x[0]);
			if (dA >= 0){
				//  not a maximum??
				return x[1];
			}
			//  the maximum is the zero of the first derivative
			dXmax = -dB / (2 * dA);
			//  check the result
			if (dXmax > x[2]){
				dXmax = x[1];
			}
			if (dXmax < x[0]){
				dXmax = x[1];
			}
			return dXmax;
		}
		/// convert coordinates from cartesian to polar
		static void c2r(RCoords& c1, CCoords& c2){
			c2r(c1, c2.x, c2.y, c2.z);
		}
		static void c2r(RCoords& c1, float x, float y, float z){
			c1.rho = sqrtf(x * x + y * y + z * z);
			c1.fi = atan2(y, x);
			c1.theta = atan2(z, sqrtf(x * x + y * y));
		}
		static float approx(float* arrx, float* arry, int range, float* coeff, int points){
			float dMa[20][20];
			float *a[20];
			float b[20];
			int i, j, k;
			float sigma, y, delta;
			//   we have memory limitations
			if (range >= 20){
				return -1.f;
			}
			//   we need more points than range+1
			if (points <= range){
				return -1.f;
			}
			//  we need more than 1 point to avoid divde by zero
			if (points <= 1){
				return -1.f;
			}
			//   fill pointers array
			for (i = 0; i < 20; i++){
				a[i] = &dMa[i][0];
				b[i] = 0.f;
			}
			//   Fill the matrix with values for polynomial approximation
			for (i = 0; i <= range; i++){
				for (j = 0; j <= range; j++){
					a[i][j] = 0.f;
					for (k = 0; k < points; k++){
						a[i][j] = a[i][j] + ipow(arrx[k], (i + j));
					}
				}
				b[i] = 0.f;
				for (k = 0; k < points; k++)
				{
					b[i] = b[i] + arry[k] * ipow(arrx[k], i);
				}
			}

			//   solve the system of linear equations
			syst_gaus(range + 1, a, b, coeff);
			//   calculate the standard deviation of the approximation
			sigma = 0.f;
			for (i = 0; i < points; i++){
				//   use Horner algorithm for polinomial calculation
				y = coeff[range];
				for (j = range; j > 0; j--){
					y = y * arrx[i] + coeff[j - 1];
				}
				delta = y - arry[i];
				sigma += delta * delta;
			}
			sigma = sqrtf(sigma / (points - 1));
			//   leave
			return sigma;
		}

		static float ipow(float x, int i){
			float p, x1;
			int j;
			if (i < 0){
				i = -i;
				x1 = 1.f / x;
			}
			else{
				x1 = x;
			}
			p = 1.f;
			for (j = 0; j < i; j++){
				p = p * x1;
			}
			return p;
		}

		static void syst_gaus(int mm, float **matr_x, float *vect_y, float *param){
			float x;
			int i, j, k;
			//---------------------------------------------
			//  Calculate the triangular matrix - matr_x.
			//---------------------------------------------
			for (i = 0; i < mm; i++){
				for (k = i + 1; k < mm; k++){
					x = matr_x[k][i] / matr_x[i][i];
					for (j = i + 1; j < mm; j++){
						matr_x[k][j] = matr_x[k][j] - matr_x[i][j] * x;
					}
					vect_y[k] = vect_y[k] - vect_y[i] * x;
				}
			}
			//-----------------------------------
			//  Back moving - Jordan method.
			//-----------------------------------
			for (i = mm - 1; i > 0; i--){
				for (k = i - 1; k >= 0; k--){
					x = matr_x[k][i] / matr_x[i][i];
					vect_y[k] = vect_y[k] - vect_y[i] * x;
				}
			}
			//-----------------------------------
			//  Calculate unknown parameters.
			//-----------------------------------
			for (i = 0; i < mm; i++){
				param[i] = vect_y[i] / matr_x[i][i];
			}
		}
	};
}



#endif /* UTILS_H_ */
