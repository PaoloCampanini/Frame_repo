// Frame.cpp
// Author: Paolo Campanini

#include "Frame.h"

#include <cmath>
#include <iostream>

#define co std::cout
#define el std::endl

Frame::Point::Point() // TESTED
{
	this->x = 0.0;
	this->y = 0.0;
	this->z = 0.0;
}

Frame::Point::Point(double _x, double _y, double _z) // TESTED
{
	this->x = _x;
	this->y = _y;
	this->z = _z;
}

Frame::Point_S::Point_S() // TESTED
{
	this->azimuth = 0.0;
	this->zenith = 0.0;
	this->distance = 0.0;
}

Frame::Point_S::Point_S(double _azimuth, double _zenith, double _distance) // TESTED
{
	this->azimuth = _azimuth;
	this->zenith = _zenith;
	this->distance = _distance;
}

Frame::Vector::Vector() // TESTED
{
	this->Vx = 0.0;
	this->Vy = 0.0;
	this->Vz = 0.0;
}

Frame::Vector::Vector(double _Vx, double _Vy, double _Vz) // TESTED
{
	this->Vx = _Vx;
	this->Vy = _Vy;
	this->Vz = _Vz;
}

Frame::R::R() // TESTED
{
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			if (i == j)
				this->r[i][j] = 1.0;
			else
				this->r[i][j] = 0.0;
		}
	}
}

Frame::A::A() // TESTED
{
	Frame::Point _P(0.0, 0.0, 0.0);
	Frame::R _R = Frame::nullRM();

	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			this->a[i][j] = _R.r[i][j];
		}
	}

	this->a[0][3] = _P.x;
	this->a[1][3] = _P.y;
	this->a[2][3] = _P.z;

	this->a[3][0] = 0.0;
	this->a[3][1] = 0.0;
	this->a[3][2] = 0.0;
	this->a[3][3] = 1.0;
}

double Frame::pi() // TESTED
{
	return 3.1415926535897931;
}

double Frame::pi2() // TESTED
{
	return Frame::pi() / 2.0;
}

double Frame::cos(double _rad)
{
	return std::cos(_rad);
}

double Frame::sin(double _rad)
{
	return std::sin(_rad);
}

double Frame::tan(double _rad)
{
	return std::tan(_rad);
}

double Frame::acos(double _x)
{
	return std::acos(_x);
}

double Frame::asin(double _x)
{
	return std::asin(_x);
}

double Frame::atan(double _x)
{
	return std::atan(_x);
}

double Frame::atan2(double _y, double _x)
{
	return std::atan2(_y, _x);
}

double Frame::pow(double _base, double _exp)
{
	return std::pow(_base, _exp);
}

double Frame::pow2(double _x)
{
	return std::pow(_x, 2);
}

double Frame::sqrt(double _num)
{
	return std::sqrt(_num);
}

double Frame::to_deg(double _rad)
{
	return (180.0 / Frame::pi()) * _rad;
}

double Frame::to_rad(double _deg)
{
	return (Frame::pi() / 180.0) * _deg;
}

Frame::Point_S Frame::to_Spherical(Point* _P) // TESTED
{
	double _distance = Frame::sqrt(Frame::pow2(_P->x) + Frame::pow2(_P->y) + Frame::pow2(_P->z));
	if (_distance == 0.0)
		return Frame::Point_S(0.0, 0.0, 0.0); // Non calcolo azimuth e zenith

	return Frame::Point_S(Frame::atan2(_P->y, _P->x), Frame::acos(_P->z / _distance), _distance);
}

Frame::Point Frame::to_Cartesian(Point_S* _Ps) // TESTED
{
	return Frame::Point(_Ps->distance * Frame::sin(_Ps->zenith) * Frame::cos(_Ps->azimuth), _Ps->distance * Frame::sin(_Ps->zenith) * Frame::sin(_Ps->azimuth), _Ps->distance * Frame::cos(_Ps->zenith));
}

Frame::Vector Frame::to_Vector(Point* _P) // TESTED
{
	return Frame::Vector(_P->x, _P->y, _P->z);
}

Frame::Point Frame::to_Point(Vector* _V) // TESTED
{
	return Frame::Point(_V->Vx, _V->Vy, _V->Vz);
}

Frame::R Frame::nullRM() // TESTED
{
	Frame::R _R;
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			if (i == j)
				_R.r[i][j] = 1.0;
			else
				_R.r[i][j] = 0.0;
		}
	}

	return _R;
}

Frame::R Frame::transpose(R _R) // TESTED
{
	Frame::R _RT;
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			_RT.r[i][j] = _R.r[j][i];
		}
	}

	return _RT;
}

Frame::A Frame::composeHTM(R _R, Point _P) // TESTED
{
	Frame::A _A;

	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			_A.a[i][j] = _R.r[i][j];
		}
	}

	_A.a[0][3] = _P.x;
	_A.a[1][3] = _P.y;
	_A.a[2][3] = _P.z;

	_A.a[3][0] = 0.0;
	_A.a[3][1] = 0.0;
	_A.a[3][2] = 0.0;
	_A.a[3][3] = 1.0;

	return _A;
}

Frame::Vector Frame::add(Vector* _V1, Vector* _V2) // TESTED
{
	return Frame::Vector(_V1->Vx + _V2->Vx, _V1->Vy + _V2->Vy, _V1->Vz + _V2->Vz);
}

Frame::Vector Frame::substruct(Vector* _V1, Vector* _V2) // TESTED
{
	return Frame::Vector(_V1->Vx - _V2->Vx, _V1->Vy - _V2->Vy, _V1->Vz - _V2->Vz);
}

double Frame::scalProd(Vector* _V1, Vector* _V2) // TESTED
{
	return _V1->Vx * _V2->Vx + _V1->Vy * _V2->Vy + _V1->Vz * _V2->Vz;	
}

Frame::Vector Frame::crossProd(Vector* _V1, Vector* _V2) // TESTED
{
	return Frame::Vector(_V1->Vy * _V2->Vz - _V1->Vz * _V2->Vy, _V1->Vz * _V2->Vx - _V1->Vx * _V2->Vz, _V1->Vx * _V2->Vy - _V1->Vy * _V2->Vx);
}

double Frame::norm(Vector* _V) // TESTED
{
	return Frame::sqrt(Frame::scalProd(_V, _V));
}

Frame::Vector Frame::normalize(Vector* _V) // TESTED
{
	double _norm = Frame::norm(_V);
	if (_norm == 0.0) // Se la norma è zero allora il vettore è un vettore nullo e ritorno un vettore nullo
	{
		return *_V;
	}
	return Frame::Vector(_V->Vx / _norm, _V->Vy / _norm, _V->Vz / _norm);
}

double Frame::cos(Vector* _V1, Vector* _V2) // TESTED
{
	double _normProd = Frame::norm(_V1) * Frame::norm(_V2);
	if (_normProd == 0.0) // Se il prodotto delle due norme è zero allora almeno uno dei due vettori è un vettore nullo e non posso calcolare l'angolo
	{
		return 1.0; // perchè acos(1.0) = 0.0
	}
	return Frame::scalProd(_V1, _V2) / _normProd;
}

double Frame::angle(Vector* _V1, Vector* _V2) // TESTED
{	
	return Frame::acos(Frame::cos(_V1, _V2)); 
}

double Frame::distance(Point* _P1, Point* _P2) // TESTED
{
	return Frame::sqrt(Frame::pow2(_P2->x - _P1->x) + Frame::pow2(_P2->y - _P1->y) + Frame::pow2(_P2->z - _P1->z));
}

double Frame::distance(Vector* _V1, Vector* _V2) // TESTED
{
	return Frame::sqrt(Frame::pow2(_V2->Vx - _V1->Vx) + Frame::pow2(_V2->Vy - _V1->Vy) + Frame::pow2(_V2->Vz - _V1->Vz));
}

Frame::Vector Frame::fromPoints(Point* _P1, Point* _P2) // TESTED
{
	return Frame::Vector(_P2->x - _P1->x, _P2->y - _P1->y, _P2->z - _P1->z);
}

void Frame::print(Point _P, disUnit _unit) // TESTED
{
	if (_unit == disUnit::noUnit)
	{
		co << "X = " << _P.x << " Y = " << _P.y << " Z = " << _P.z << el;
		return;
	}

	if (_unit == disUnit::m)
	{
		co << "X = " << _P.x << " m  Y = " << _P.y << " m  Z = " << _P.z << " m" << el;
		return;
	}

	if (_unit == disUnit::mm)
	{
		co << "X = " << _P.x * 1000.0 << " mm  Y = " << _P.y * 1000.0 << " mm  Z = " << _P.z * 1000.0 << " mm" << el;
		return;
	}

	if (_unit == disUnit::micron)
	{
		co << "X = " << _P.x * 1000000.0 << " micron  Y = " << _P.y * 1000000.0 << " micron  Z = " << _P.z * 1000000.0 << " micron" << el;
		return;
	}
}

void Frame::print(Point_S _Ps, bool _deg) // TESTED
{
	if(_deg)
		co << "Azimuth = " << Frame::to_deg(_Ps.azimuth) << " deg  Zenith = " << Frame::to_deg(_Ps.zenith) << " deg  Distance = " << _Ps.distance << el;
	else
		co << "Azimuth = " << _Ps.azimuth << " rad  Zenith = " << _Ps.zenith << " rad  Distance = " << _Ps.distance << el;
}

void Frame::print(Vector _V, disUnit _unit) // TESTED
{
	if (_unit == disUnit::noUnit)
	{
		co << "Vx = " << _V.Vx << " Vy = " << _V.Vy << " Vz = " << _V.Vz << el;
		return;
	}

	if (_unit == disUnit::m)
	{
		co << "Vx = " << _V.Vx << " m  Vy = " << _V.Vy << " m  Vz = " << _V.Vz << " m" << el;
		return;
	}

	if (_unit == disUnit::mm)
	{
		co << "Vx = " << _V.Vx * 1000.0 << " mm  Vy = " << _V.Vy * 1000.0 << " mm  Vz = " << _V.Vz * 1000.0 << " mm" << el;
		return;
	}

	if (_unit == disUnit::micron)
	{
		co << "Vx = " << _V.Vx * 1000000.0 << " micron  Vy = " << _V.Vy * 1000000.0 << " micron  Vz = " << _V.Vz * 1000000.0 << " micron" << el;
		return;
	}
}

void Frame::print(R _R) // TESTED
{
	for (int j = 0; j < 3; j++)
	{
		for (int i = 0; i < 3; i++) // Prima stampo le righe e poi le colonne
		{
			co << _R.r[i][j] << " ";
		}
		co << el;
	}
}

void Frame::print(A _A) // TESTED
{
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			co << _A.a[i][j] << " ";
		}
		co << el;
	}
}