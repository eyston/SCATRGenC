#ifndef DRIFT_H_
#define DRIFT_H_

#include "structures.h"

void drift_tp(const size_t ntp, const float mass_sun, vec3_scalar_t &h, const vec3_scalar_t &vh, const float dt, int **istat);

#endif /* DRIFT_H_ */