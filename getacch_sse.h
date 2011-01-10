#ifndef GETACHH_SSE_H_
#define GETACHH_SSE_H_

#include <xmmintrin.h>
#include "structures.h"

void getacch_sse(const size_t nbod, __m128 mass[NPLMAX], const __m128 xj[NPLMAX], const __m128 yj[NPLMAX], const __m128 zj[NPLMAX], 
	const __m128 xh[NPLMAX], const __m128 yh[NPLMAX], const __m128 zh[NPLMAX], __m128 axh[NPLMAX], __m128 ayh[NPLMAX], __m128 azh[NPLMAX]);

void xgetacch_sse(const size_t nbod, planets_sse_t & __restrict p);

void x_getacch_ir3_sse_test(const size_t nbod, const vec3_sse_t &pos, vec1_sse_t &ir3);

void x_getacch_ah0_sse_test(const size_t nbod, const vec1_sse_t mass, const vec1_sse_t ir3h, const vec3_sse_t &h, float &axh0f, float &ayh0f, float &azh0f);
void x_getacch_ah1_sse_test(const size_t nbod, const vec1_sse_t &mass, const vec1_sse_t &ir3h, const vec1_sse_t &ir3j, const vec3_sse_t &h, const vec3_sse_t &j, vec3_sse_t &ah1);
void x_getacch_ah2_sse_test(const size_t nbod, const vec1_sse_t &mass, const vec1_sse_t &ir3j, const vec3_sse_t &j, vec3_sse_t &ah2);
void x_getacch_ah3_sse_test(const size_t nbod, vec1_sse_t &mass, const vec3_sse_t &h, vec3_sse_t &ah3);

#endif /* GETACHH_SSE_H_ */
