/*
 * Vector3R.h
 *
 * obs.: inspired on the corresponding EvtGen class.
 *
 * Created on : Feb 25, 2016
 *      Author: Antonio Augusto Alves Junior
 */

/*
 *   This file is part of MCBooster.
 *
 *   MCBooster is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   MCBooster is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with MCBooster.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef VECTOR3R_H_
#define VECTOR3R_H_



#include <mcbooster/Config.h>
#include <mcbooster/GTypes.h>
#include <iosfwd>
#include <iostream>
#include <math.h>

using std::ostream;
namespace MCBooster
{
class Vector3R
{

	__host__      __device__      friend Vector3R rotateEuler(const Vector3R& v,
			GReal_t phi, GReal_t theta, GReal_t ksi);

	__host__      __device__      inline friend Vector3R operator*(GReal_t c,
			const Vector3R& v2);
	__host__      __device__      inline friend GReal_t operator*(const Vector3R& v1,
			const Vector3R& v2);
	__host__      __device__      inline friend Vector3R operator+(const Vector3R& v1,
			const Vector3R& v2);
	__host__      __device__      inline friend Vector3R operator-(const Vector3R& v1,
			const Vector3R& v2);
	__host__      __device__      inline friend Vector3R operator*(const Vector3R& v1,
			GReal_t c);
	__host__      __device__      inline friend Vector3R operator/(const Vector3R& v1,
			GReal_t c);
	__host__      __device__      friend Vector3R cross(const Vector3R& v1,
			const Vector3R& v2);

public:
	__host__ __device__ inline Vector3R();
	__host__ __device__ inline Vector3R(GReal_t x, GReal_t y, GReal_t z);
	__host__ __device__ inline Vector3R(const Vector3R& other);
	__host__      __device__      inline Vector3R& operator*=(const GReal_t c);
	__host__      __device__      inline Vector3R& operator/=(const GReal_t c);
	__host__      __device__      inline Vector3R& operator+=(const Vector3R& v2);
	__host__      __device__      inline Vector3R& operator-=(const Vector3R& v2);
	__host__ __device__ inline void set(GInt_t i, GReal_t d);
	__host__ __device__ inline void set(GReal_t x, GReal_t y, GReal_t z);
	__host__ __device__ inline void applyRotateEuler(GReal_t phi, GReal_t theta,
			GReal_t ksi);
	__host__      __device__      inline GReal_t get(GInt_t i) const;
	__host__       inline friend std::ostream& operator<<(std::ostream& s,
			const Vector3R& v);
	__host__      __device__      inline GReal_t dot(const Vector3R& v2);
	__host__      __device__      inline GReal_t d3mag() const;

private:

	GReal_t v[3];

};

inline Vector3R& Vector3R::operator*=(const GReal_t c)
{

	v[0] *= c;
	v[1] *= c;
	v[2] *= c;
	return *this;
}

inline Vector3R& Vector3R::operator/=(const GReal_t c)
{

	v[0] /= c;
	v[1] /= c;
	v[2] /= c;
	return *this;
}

inline Vector3R& Vector3R::operator+=(const Vector3R& v2)
{

	v[0] += v2.v[0];
	v[1] += v2.v[1];
	v[2] += v2.v[2];
	return *this;
}

inline Vector3R& Vector3R::operator-=(const Vector3R& v2)
{

	v[0] -= v2.v[0];
	v[1] -= v2.v[1];
	v[2] -= v2.v[2];
	return *this;
}

inline Vector3R operator*(GReal_t c, const Vector3R& v2)
{

	return Vector3R(v2) *= c;
}

inline Vector3R operator*(const Vector3R& v1, GReal_t c)
{

	return Vector3R(v1) *= c;
}

inline Vector3R operator/(const Vector3R& v1, GReal_t c)
{

	return Vector3R(v1) /= c;
}

inline GReal_t operator*(const Vector3R& v1, const Vector3R& v2)
{

	return v1.v[0] * v2.v[0] + v1.v[1] * v2.v[1] + v1.v[2] * v2.v[2];
}

inline Vector3R operator+(const Vector3R& v1, const Vector3R& v2)
{

	return Vector3R(v1) += v2;
}

inline Vector3R operator-(const Vector3R& v1, const Vector3R& v2)
{

	return Vector3R(v1) -= v2;

}

inline GReal_t Vector3R::get(GInt_t i) const
{
	return v[i];
}

inline void Vector3R::set(GInt_t i, GReal_t d)
{

	v[i] = d;
}

inline void Vector3R::set(GReal_t x, GReal_t y, GReal_t z)
{

	v[0] = x;
	v[1] = y;
	v[2] = z;
}

inline Vector3R::Vector3R()
{

	v[0] = v[1] = v[2] = 0.0;
}

inline Vector3R::Vector3R(GReal_t x, GReal_t y, GReal_t z)
{

	v[0] = x;
	v[1] = y;
	v[2] = z;
}

inline Vector3R::Vector3R(const Vector3R& other)
{

	v[0] = other.get(0);
	v[1] = other.get(1);
	v[2] = other.get(2);
}

inline Vector3R rotateEuler(const Vector3R& v, GReal_t alpha, GReal_t beta,
		GReal_t gamma)
{

	Vector3R tmp(v);
	tmp.applyRotateEuler(alpha, beta, gamma);
	return tmp;

}

inline void Vector3R::applyRotateEuler(GReal_t phi, GReal_t theta, GReal_t ksi)
{

	GReal_t temp[3];
	GReal_t sp, st, sk, cp, ct, ck;

	sp = sin(phi);
	st = sin(theta);
	sk = sin(ksi);
	cp = cos(phi);
	ct = cos(theta);
	ck = cos(ksi);

	temp[0] = (ck * ct * cp - sk * sp) * v[0] + (-sk * ct * cp - ck * sp) * v[1]
			+ st * cp * v[2];
	temp[1] = (ck * ct * sp + sk * cp) * v[0] + (-sk * ct * sp + ck * cp) * v[1]
			+ st * sp * v[2];
	temp[2] = -ck * st * v[0] + sk * st * v[1] + ct * v[2];

	v[0] = temp[0];
	v[1] = temp[1];
	v[2] = temp[2];
}

inline ostream& operator<<(ostream& s, const Vector3R& v)
{

	s << "(" << v.v[0] << "," << v.v[1] << "," << v.v[2] << ")";

	return s;

}

inline Vector3R cross(const Vector3R& p1, const Vector3R& p2)
{

	//Calcs the cross product.  Added by djl on July 27, 1995.
	//Modified for real vectros by ryd Aug 28-96

	return Vector3R(p1.v[1] * p2.v[2] - p1.v[2] * p2.v[1],
			p1.v[2] * p2.v[0] - p1.v[0] * p2.v[2],
			p1.v[0] * p2.v[1] - p1.v[1] * p2.v[0]);

}

inline GReal_t Vector3R::d3mag() const

// returns the 3 momentum mag.
{
	GReal_t temp;

	temp = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
	temp = sqrt(temp);

	return temp;
} // r3mag

inline GReal_t Vector3R::dot(const Vector3R& p2)
{

	GReal_t temp;

	temp = v[0] * p2.v[0];
	temp += v[0] * p2.v[0];
	temp += v[0] * p2.v[0];

	return temp;
} //dot
}
#endif /* VECTOR3R_H_ */
