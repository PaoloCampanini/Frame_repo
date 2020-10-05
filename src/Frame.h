// Frame.h
// Author: Paolo Campanini

#ifndef _FRAME_H
#define _FRAME_H

namespace Frame
{
	enum disUnit
	{
		noUnit, // no unit
		micron, // micrometri (1e-6 m)
		mm,	// millimetri (1e-3 m)
		m // metri
	}; // Units of measurement for distances

	class Point
	{
	public:
		double x; // [m]
		double y; // [m]
		double z; // [m]

		Point(); // Origin (0, 0, 0)
		Point(double _x, double _y, double _z);
	}; // Cartesian coordinates of a point in space

	class Point_S
	{
	public:
		double azimuth;
		double zenith;
		double distance;

		Point_S(); // Origin (0, 0, 0)
		Point_S(double _azimuth, double _zenith, double _distance);
	}; // Spherical coordinates of a point in space

	class Vector
	{
	public:
		double Vx; // X component
		double Vy; // Y component
		double Vz; // Z component

		Vector(); // Null vector [0, 0, 0]
		Vector(double _Vx, double _Vy, double _Vz);
	};

	class R
	{
	public:
		double r[3][3];

		R(); // Null rotation (identity matrix)
	}; // Rotation matrix

	class A
	{
	public:
		double a[4][4];

		A(); // Composition of null rotation (identity matrix) and origin (0, 0, 0)
	}; // Homogeneous transformation matrix

// Methods
	// Trigonometric functions
	double cos(double _rad); // (argument in radiants and must be in [-1,+1]).
	double sin(double _rad); // (argument in radiants and must be in [-1,+1]).
	double tan(double _rad); // (argument in radiants).
	double acos(double _x); // (return an angle in radiants). (input in interval [-1, 1], output in interval [0, +pi]).
	double asin(double _x); // (return an angle in radiants). (input in interval [-1, 1], output in interval [-pi/2, +pi/2]).
	double atan(double _x); // (return an angle in radiants in [-pi/2,+pi/2]).
	double atan2(double _y, double _x); // (return an angle in radiants in [-pi,+pi]).

	// Power functions
	double pow(double _base, double _exp);
	double pow2(double _x); // (return _x^2).
	double sqrt(double _num);

	// Matrix
	double pi(); // Return 'pi' (3.141...).
	double pi2(); // Return 'pi/2'.

	//* (pointers are sometime used to reduce the amount of memory when passing parameters to functions)
	// Conversions
	double to_deg(double _rad); // Convert from radian to degree.
	double to_rad(double _deg); // Convert from degree to radian.
	Point_S to_Spherical(Point* _P); // Conversion from Cartesian to Spherical coordinates. (azimuth [-pi, +pi], zenith [0, +pi], distance>=0).
	Point to_Cartesian(Point_S* _Ps); // Conversion from Spherical to Cartesian coordinates.
	Vector to_Vector(Point* _P); // Converte a point object into a vector object.
	Point to_Point(Vector* _V); // Converte a vector object into a point object.

	R nullRM(); // Return a rotation matrix that does not rotate a frame.
	R transpose(R _R); // Return the transpose of a rotation matrix.
	A composeHTM(R _R, Point _P); // Return an HTM with the indicated arguments.

	Vector add(Vector* _V1, Vector* _V2); // Return the sum of two vectors (V = V1+V2).
	Vector substruct(Vector* _V1, Vector* _V2); // Return the difference of two vectors (V = V1-V2).
	double scalProd(Vector* _V1, Vector* _V2); // Return the scalar product of two vectors.
	Vector crossProd(Vector* _V1, Vector* _V2); // Return the cross product of two vectors.
	double norm(Vector* _V); // Return the 2-norm of a vector.
	Vector normalize(Vector* _V); // Return normalized vector such that its norm is =1.
	double cos(Vector* _V1, Vector* _V2); // Return the cosine of the angle between two vectors.
	double angle(Vector* _V1, Vector* _V2); // Return the angle in radian between two vectors. The angle can be in the range [0, pi].
	double distance(Point* _P1, Point* _P2); // Return the Euclidean distance between two points in space.
	double distance(Vector* _V1, Vector* _V2); // Return the Euclidean distance between two points defined by two vectors.
	Vector fromPoints(Point* _P1, Point* _P2); // Return the vector from 'P1' to 'P2'; i.e. 'P2-P1'.

	void print(Point _P, disUnit _unit = disUnit::noUnit); // Print on console a point in cartesian coordinates.
	void print(Point_S _Ps, bool _deg = true); // Print on console a point in spherical coordinates. Argument '_deg' determines the unit of the angles.
	void print(Vector _V, disUnit _unit = disUnit::noUnit); // Print on console the components of a vector.
	void print(R _R); // Print on console a rotation matrix.
	void print(A _A); // Print on console an homogeneous transformation matrix.

} // namespace Frame

#endif //_FRAME_H