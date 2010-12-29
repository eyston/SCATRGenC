#ifndef GETACHH_SSE_H_
#define GETACHH_SSE_H_

#include <xmmintrin.h>
#include "structures.h"

void getacch_sse(const size_t nbod, __m128 mass[NPLMAX], const __m128 xj[NPLMAX], const __m128 yj[NPLMAX], const __m128 zj[NPLMAX], 
	const __m128 xh[NPLMAX], const __m128 yh[NPLMAX], const __m128 zh[NPLMAX], __m128 axh[NPLMAX], __m128 ayh[NPLMAX], __m128 azh[NPLMAX]);

void xgetacch_sse(const size_t nbod, planets_sse_t & __restrict p);


#endif /* GETACHH_SSE_H_ */
