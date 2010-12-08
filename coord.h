/*
 * coord.h
 *
 *  Created on: Nov 23, 2010
 *      Author: Huey
 */

#ifndef COORD_H_
#define COORD_H_

#include "structures.h"

void coord_h2j(NBodies bodies);

void coord_j2h(NBodies bodies);

void coord_h2j(const size_t nbod, const float mass[NPLMAX], const float xh[NPLMAX], const float yh[NPLMAX], const float zh[NPLMAX], const float vxh[NPLMAX], const float vyh[NPLMAX], const float vzh[NPLMAX],
	float xj[NPLMAX], float yj[NPLMAX], float zj[NPLMAX], float vxj[NPLMAX], float vyj[NPLMAX], float vzj[NPLMAX]);

#endif /* COORD_H_ */
