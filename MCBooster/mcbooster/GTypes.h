/*
 * GTypes.h
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


#ifndef GTYPES_H_
#define GTYPES_H_

namespace MCBooster
{
//---- types -------------------------------------------------------------------

typedef char GChar_t;      ///< Signed Character 1 byte (char)
typedef unsigned char GUChar_t;  ///< Unsigned Character 1 byte (unsigned char)
typedef short GShort_t;     ///< Signed Short integer 2 bytes (short)
typedef unsigned short GUShort_t; ///<Unsigned Short integer 2 bytes (unsigned short)
typedef int GInt_t;       ///< Signed integer 4 bytes (int)
typedef unsigned int GUInt_t;      ///< Unsigned integer 4 bytes (unsigned int)
typedef long GLong_t;      ///< Signed long integer 4 bytes (long)
typedef unsigned long GULong_t; //Unsigned long integer 4 bytes (unsigned long)
typedef float GFloat_t;     ///< Float 4 bytes (float)
typedef double GDouble_t;    ///< Double 8 bytes
typedef long double GLongDouble_t;    ///< Long Double
typedef char GText_t;      ///< General string (char)
typedef bool GBool_t;      ///< Boolean (0=false, 1=true) (bool)
typedef unsigned char GByte_t;      ///< Byte (8 bits) (unsigned char)
typedef long long GLong64_t; ///< Portable signed long integer 8 bytes
typedef unsigned long long GULong64_t; ///< Portable unsigned long integer 8 bytes
#ifdef FP_SINGLE
typedef float GReal_t;///< Double 8 bytes or float 4 bytes
#else
typedef double GReal_t;///< Double 8 bytes or float 4 bytes
#endif

//---- constants ---------------------------------------------------------------

#ifndef NULL
#define NULL 0
#endif

const GBool_t kTrue = true;
const GBool_t kFalse = false;


#define kMAXP  9
#define PI     3.1415926535897932384626422832795028841971

}
#endif /* GTYPES_H_ */
