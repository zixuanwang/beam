#ifndef COORDS_H_
#define COORDS_H_

namespace Beam{
	/// Polar coordinate system.
	class RCoords {
	public:
		float fi;
		float theta;
		float rho;
	};
	/// Cartesian coordinate system.
	class CCoords{
	public:
		float x;
		float y;
		float z;
	};
}



#endif /* COORDS_H_ */
