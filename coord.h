#ifndef COORD_H_
#define COORD_H_

#include "structures.h"

void coord_h2j(const size_t nbod, const float mass[NPLMAX], const float xh[NPLMAX], const float yh[NPLMAX], const float zh[NPLMAX], const float vxh[NPLMAX], const float vyh[NPLMAX], const float vzh[NPLMAX],
	float xj[NPLMAX], float yj[NPLMAX], float zj[NPLMAX], float vxj[NPLMAX], float vyj[NPLMAX], float vzj[NPLMAX]);

void xcoord_h2j(const size_t nbod, planets_scalar_t &planets);
#endif /* COORD_H_ */
