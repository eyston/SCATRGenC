#ifndef GETACHH_H_
#define GETACHH_H_

#include <xmmintrin.h>

void test_getacch_sse(const size_t nbod, size_t nbod_v, __m128 mass[NPLMAX], const __m128 xj[NPLMAX], const __m128 yj[NPLMAX], const __m128 zj[NPLMAX], const __m128 xh[NPLMAX], const __m128 yh[NPLMAX], const __m128 zh[NPLMAX]);

void getacch(const size_t nbod, const float mass[NPLMAX], const float xj[NPLMAX], const float yj[NPLMAX], const float zj[NPLMAX], 
	const float xh[NPLMAX], const float yh[NPLMAX], const float zh[NPLMAX], float axh[NPLMAX], float ayh[NPLMAX], float azh[NPLMAX]);

static void getacch_ir3(const size_t nbod, const size_t istart, const float x[NPLMAX], const float y[NPLMAX], const float z[NPLMAX], float ir3[NPLMAX]);

static void getacch_ah0(const size_t istart, const size_t nbod, const float mass[NPLMAX], const float xh[NPLMAX], const float yh[NPLMAX], const float zh[NPLMAX], 
	const float ir3h[NPLMAX], float &axh0, float &ayh0, float &azh0);

static void getacch_ah1(const size_t nbod, const float mass[NPLMAX], const float xh[NPLMAX], const float yh[NPLMAX], const float zh[NPLMAX], const float xj[NPLMAX], const float yj[NPLMAX], const float zj[NPLMAX], 
	const float ir3h[NPLMAX], const float ir3j[NPLMAX], float axh1[NPLMAX], float ayh1[NPLMAX], float azh1[NPLMAX]);

static void getacch_ah2(const size_t nbod, const float mass[NPLMAX], const float xj[NPLMAX], const float yj[NPLMAX], const float zj[NPLMAX], const float ir3j[NPLMAX], 
	float axh2[NPLMAX], float ayh2[NPLMAX], float azh2[NPLMAX]);

static void getacch_ah3(const size_t nbod, const float mass[NPLMAX], const float xh[NPLMAX], const float yh[NPLMAX], const float zh[NPLMAX], float axh3[NPLMAX], float ayh3[NPLMAX], float azh3[NPLMAX]);

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


#endif /* GETACHH_H_ */
