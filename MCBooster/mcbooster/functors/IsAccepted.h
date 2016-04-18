
/*
 * IsAccepted.h
 *
 * Copyright 2016 Antonio Augusto Alves Junior
 *
 * Created on : 29/03/2016
 *      Author: Antonio Augusto Alves Junior
 */

/*
    This file is part of MCBooster.

    MCBooster is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    MCBooster is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with MCBooster.  If not, see <http://www.gnu.org/licenses/>.
*/

/**\file IsAccepted.h
 * Implements isAccepted.
 */

#ifndef ISACCEPTED_H_
#define ISACCEPTED_H_


#include <mcbooster/Config.h>
#include <mcbooster/GTypes.h>

namespace MCBooster
{

struct isAccepted
{
  __host__ __device__
  bool operator()(const int x)
  {
    return (x == 1 ) ;
  }
};


}


#endif /* ISACCEPTED_H_ */
