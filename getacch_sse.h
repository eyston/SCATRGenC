#ifndef GETACHH_SSE_H_
#define GETACHH_SSE_H_

#include <xmmintrin.h>

void getacch_sse(const size_t nbod, __m128 mass[NPLMAX], const __m128 xj[NPLMAX], const __m128 yj[NPLMAX], const __m128 zj[NPLMAX], 
	const __m128 xh[NPLMAX], const __m128 yh[NPLMAX], const __m128 zh[NPLMAX], __m128 axh[NPLMAX], __m128 ayh[NPLMAX], __m128 azh[NPLMAX]);

static void getacch_ir3_sse(const size_t nbod, const __m128 x[NPLMAX], const __m128 y[NPLMAX], const __m128 z[NPLMAX], __m128 ir3[NPLMAX]);

static void getacch_ah0_sse(const size_t nbod, const __m128 mass[NPLMAX], const __m128 xh[NPLMAX], const __m128 yh[NPLMAX], const __m128 zh[NPLMAX], 
	const __m128 ir3h[NPLMAX], __m128 &axh0, __m128 &ayh0, __m128 &azh0);

static void getacch_ah1_sse(const size_t nbod, const __m128 mass[NPLMAX], const __m128 xh[NPLMAX], const __m128 yh[NPLMAX], const __m128 zh[NPLMAX], const __m128 xj[NPLMAX], const __m128 yj[NPLMAX], const __m128 zj[NPLMAX], 
	const __m128 ir3h[NPLMAX], const __m128 ir3j[NPLMAX], __m128 axh1[NPLMAX], __m128 ayh1[NPLMAX], __m128 azh1[NPLMAX]);

static void getacch_ah2_sse(const size_t nbod, const __m128 mass[NPLMAX], const __m128 xj[NPLMAX], const __m128 yj[NPLMAX], const __m128 zj[NPLMAX], const __m128 ir3j[NPLMAX], 
	__m128 axh2[NPLMAX], __m128 ayh2[NPLMAX], __m128 azh2[NPLMAX]);

static void getacch_ah3_sse(const size_t nbod, __m128 mass[NPLMAX], const __m128 xh[NPLMAX], const __m128 yh[NPLMAX], const __m128 zh[NPLMAX], __m128 axh3[NPLMAX], __m128 ayh3[NPLMAX], __m128 azh3[NPLMAX]);

#endif /* GETACHH_SSE_H_ */
