#include <stdio.h>
#include <xmmintrin.h>
#include <math.h>
#include "structures.h"
#include "getacch.h"

static float *ir3h_s = (float*) _mm_malloc(NPLMAX * sizeof(float), 16);
static float *irh_s = (float*) _mm_malloc(NPLMAX * sizeof(float), 16);
static float *ir3j_s = (float*) _mm_malloc(NPLMAX * sizeof(float), 16);
static float *irj_s = (float*) _mm_malloc(NPLMAX * sizeof(float), 16);
static __m128 *ir3h_v = (__m128*) ir3h_s;
static __m128 *irh_v = (__m128*) irh_s;
static __m128 *ir3j_v = (__m128*) ir3j_s;
static __m128 *irj_v = (__m128*) irj_s;

static float *axh1_s = (float*) _mm_malloc(NPLMAX * sizeof(float), 16);
static float *ayh1_s = (float*) _mm_malloc(NPLMAX * sizeof(float), 16);
static float *azh1_s = (float*) _mm_malloc(NPLMAX * sizeof(float), 16);
static __m128 *axh1_v = (__m128*) axh1_s;
static __m128 *ayh1_v = (__m128*) ayh1_s;
static __m128 *azh1_v = (__m128*) azh1_s;

static float *axh2_s = (float*) _mm_malloc(NPLMAX * sizeof(float), 16);
static float *ayh2_s = (float*) _mm_malloc(NPLMAX * sizeof(float), 16);
static float *azh2_s = (float*) _mm_malloc(NPLMAX * sizeof(float), 16);
static __m128 *axh2_v = (__m128*) axh2_s;
static __m128 *ayh2_v = (__m128*) ayh2_s;
static __m128 *azh2_v = (__m128*) azh2_s;

static float *axh3_s = (float*) _mm_malloc(NPLMAX * sizeof(float), 16);
static float *ayh3_s = (float*) _mm_malloc(NPLMAX * sizeof(float), 16);
static float *azh3_s = (float*) _mm_malloc(NPLMAX * sizeof(float), 16);
static __m128 *axh3_v = (__m128*) axh3_s;
static __m128 *ayh3_v = (__m128*) ayh3_s;
static __m128 *azh3_v = (__m128*) azh3_s;

void test_getacch_sse(const size_t nbod, const size_t nbod_v, __m128 mass[NPLMAX], const __m128 xj[NPLMAX], const __m128 yj[NPLMAX], const __m128 zj[NPLMAX], const __m128 xh[NPLMAX], const __m128 yh[NPLMAX], const __m128 zh[NPLMAX])
{
	float *mass_s = (float*)mass;
	float *xj_s = (float*)xj;
	float *yj_s = (float*)yj;
	float *zj_s = (float*)zj;
	float *xh_s = (float*)xh;
	float *yh_s = (float*)yh;
	float *zh_s = (float*)zh;

	getacch_ir3(nbod, 1, xj_s, yj_s, zj_s, ir3j_s);
	getacch_ir3(nbod, 1, xh_s, yh_s, zh_s, ir3h_s);

	printf("SCALAR\n");
	for(size_t i = 0; i < nbod; ++i)
	{
		printf("ir3j, irj: %E, %E\n", ir3j_s[i], irj_s[i]);
		printf("ir3h, irh: %E, %E\n", ir3h_s[i], irh_s[i]);
	}

	getacch_ir3_sse(nbod_v, xh, yh, zh, ir3h_v);
	getacch_ir3_sse(nbod_v, xj, yj, zj, ir3j_v);

	printf("VECTOR\n");
	for(size_t i = 0; i < nbod; ++i)
	{
		printf("ir3j, irj: %E, %E\n", ir3j_s[i], irj_s[i]);
		printf("ir3h, irh: %E, %E\n", ir3h_s[i], irh_s[i]);
	}

	float axh0_s = 0.0f, ayh0_s = 0.0f, azh0_s = 0.0f;
	const __m128 zeroes = _mm_set1_ps(0.0f);

	__m128 axh0_v = zeroes, ayh0_v = zeroes, azh0_v = zeroes;

	getacch_ah0(2, nbod, mass_s, xh_s, yh_s, zh_s, ir3h_s, axh0_s, ayh0_s, azh0_s);
	printf("SCALAR\n");
	printf("axh0, ayh0, azh0: %E, %E, %E\n", axh0_s, ayh0_s, azh0_s);

	getacch_ah0_sse(nbod_v, mass, xh, yh, zh, ir3h_v, axh0_v, ayh0_v, azh0_v);
	printf("VECTOR\n");
	printf("axh0, ayh0, azh0: %E, %E, %E\n", ((float*)&axh0_v)[0], ((float*)&ayh0_v)[0], ((float*)&azh0_v)[0]);

	getacch_ah1(nbod, mass_s, xh_s, yh_s, zh_s, xj_s, yj_s, zj_s, ir3h_s, ir3j_s, axh1_s, ayh1_s, azh1_s);
	printf("SCALAR\n");
	for(size_t i = 0; i < nbod; ++i)
	{
		printf("axh1, ayh1, azh1: %E, %E, %E\n", axh1_s[i], ayh1_s[i], azh1_s[i]);
	}
	
	getacch_ah1_sse(nbod_v, mass, xh, yh, zh, xj, yj, zj, ir3h_v, ir3j_v, axh1_v, ayh1_v, azh1_v);
	printf("VECTOR\n");
	for(size_t i = 0; i < nbod; ++i)
	{
		printf("axh1, ayh1, azh1: %E, %E, %E\n", axh1_s[i], ayh1_s[i], azh1_s[i]);
	}

	getacch_ah2(nbod, mass_s, xj_s, yj_s, zj_s, ir3j_s, axh2_s, ayh2_s, azh2_s);
	printf("SCALAR\n");
	for(size_t i = 0; i < nbod; ++i)
	{
		printf("axh2, ayh2, azh2: %E, %E, %E\n", axh2_s[i], ayh2_s[i], azh2_s[i]);
	}

	getacch_ah2_sse(nbod_v, mass, xj, yj, zj, ir3j_v, axh2_v, ayh2_v, azh2_v);
	printf("VECTOR\n");
	for(size_t i = 0; i < nbod; ++i)
	{
		printf("axh2, ayh2, azh2: %E, %E, %E\n", axh2_s[i], ayh2_s[i], azh2_s[i]);
	}


	getacch_ah3(nbod, mass_s, xh_s, yh_s, zh_s, axh3_s, ayh3_s, azh3_s);
	printf("SCALAR\n");
	for(size_t i = 0; i < nbod; ++i)
	{
		printf("axh3, ayh3, azh3: %E, %E, %E\n", axh3_s[i], ayh3_s[i], azh3_s[i]);
	}

	printf("SCALAR\n");
	printf("P3: axh: %E\n", axh0_s + axh1_s[2] + axh2_s[2] + axh3_s[2]);

	getacch_ah3_sse(nbod_v, mass, xh, yh, zh, axh3_v, ayh3_v, azh3_v);
	printf("VECTOR\n");
	for(size_t i = 0; i < nbod; ++i)
	{
		printf("axh3, ayh3, azh3: %E, %E, %E\n", axh3_s[i], ayh3_s[i], azh3_s[i]);
	}

	printf("VECTOR\n");
	printf("P3: axh: %E\n", axh0_s + axh1_s[2] + axh2_s[2] + axh3_s[2]);

}


void getacch(const size_t nbod, const float mass[NPLMAX], const float xj[NPLMAX], const float yj[NPLMAX], const float zj[NPLMAX], 
	const float xh[NPLMAX], const float yh[NPLMAX], const float zh[NPLMAX], float axh[NPLMAX], float ayh[NPLMAX], float azh[NPLMAX])
{
	float ir3j[NPLMAX];
	float ir3h[NPLMAX];

	getacch_ir3(nbod, 1, xj, yj, zj, ir3j);
	getacch_ir3(nbod, 1, xh, yh, zh, ir3h);

	float axh0 = 0.0f;
	float ayh0 = 0.0f;
	float azh0 = 0.0f;

	getacch_ah0(2, nbod, mass, xh, yh, zh, ir3h, axh0, ayh0, azh0);

	float axh1[NPLMAX];
	float ayh1[NPLMAX];
	float azh1[NPLMAX];

	getacch_ah1(nbod, mass, xh, yh, zh, xj, yj, zj, ir3h, ir3j, axh1, ayh1, azh1);

	float axh2[NPLMAX];
	float ayh2[NPLMAX];
	float azh2[NPLMAX];

	getacch_ah2(nbod, mass, xj, yj, zj, ir3j, axh2, ayh2, azh2);

	float axh3[NPLMAX];
	float ayh3[NPLMAX];
	float azh3[NPLMAX];

	getacch_ah3(nbod, mass, xh, yh, zh, axh3, ayh3, azh3);

	axh[0] = 0.0f;
	ayh[0] = 0.0f;
	azh[0] = 0.0f;

	for(size_t i = 1; i < nbod; ++i)
	{
		axh[i] = axh0 + axh1[i] + axh2[i] + axh3[i];
		ayh[i] = ayh0 + ayh1[i] + ayh2[i] + ayh3[i];
		azh[i] = azh0 + azh1[i] + azh2[i] + azh3[i];
	}

	// TODO: j2 and j4 stuff ... not implemented
}

static void getacch_ir3(const size_t nbod, const size_t istart, const float x[NPLMAX], const float y[NPLMAX], const float z[NPLMAX], float ir3[NPLMAX])
{
	for(size_t i = istart; i < nbod; ++i)
	{
		float r2 = x[i]*x[i] + y[i]*y[i] + z[i]*z[i];
		float ir = 1.0f / sqrt(r2);
		ir3[i] = ir / r2;
	}
}

static void getacch_ah0(const size_t istart, const size_t nbod, const float mass[NPLMAX], const float xh[NPLMAX], const float yh[NPLMAX], const float zh[NPLMAX], 
	const float ir3h[NPLMAX], float &axh0, float &ayh0, float &azh0)
{
	axh0 = 0.0f;
	ayh0 = 0.0f;
	azh0 = 0.0f;

	for(size_t i = istart; i < nbod; ++i)
	{
		float fac = mass[i] * ir3h[i];

		axh0 -= fac * xh[i];
		ayh0 -= fac * yh[i];
		azh0 -= fac * zh[i];
	}
}

static void getacch_ah1(const size_t nbod, const float mass[NPLMAX], const float xh[NPLMAX], const float yh[NPLMAX], const float zh[NPLMAX], const float xj[NPLMAX], const float yj[NPLMAX], const float zj[NPLMAX], 
	const float ir3h[NPLMAX], const float ir3j[NPLMAX], float axh1[NPLMAX], float ayh1[NPLMAX], float azh1[NPLMAX])
{
	axh1[0] = 0.0f;
	ayh1[0] = 0.0f;
	azh1[0] = 0.0f;

	axh1[1] = 0.0f;
	ayh1[1] = 0.0f;
	azh1[1] = 0.0f;

	for(size_t i = 2; i < nbod; ++i)
	{
		float ah1j = xj[i] * ir3j[i];
		float ah1h = xh[i] * ir3h[i];
		axh1[i] = mass[0]*(ah1j - ah1h);

		ah1j = yj[i] * ir3j[i];
		ah1h = yh[i] * ir3h[i];
		ayh1[i] = mass[0]*(ah1j - ah1h);

		ah1j = zj[i] * ir3j[i];
		ah1h = zh[i] * ir3h[i];
		azh1[i] = mass[0]*(ah1j - ah1h);
	}
}

static void getacch_ah2(const size_t nbod, const float mass[NPLMAX], const float xj[NPLMAX], const float yj[NPLMAX], const float zj[NPLMAX], const float ir3j[NPLMAX], 
	float axh2[NPLMAX], float ayh2[NPLMAX], float azh2[NPLMAX])
{
	axh2[0] = 0.0f;
	ayh2[0] = 0.0f;
	azh2[0] = 0.0f;

	axh2[1] = 0.0f;
	ayh2[1] = 0.0f;
	azh2[1] = 0.0f;

	float etaj = mass[0];
	for(size_t i = 2; i < nbod; ++i)
	{
		etaj += mass[i-1];
		float fac = mass[i] * mass[0] * ir3j[i] / etaj;
		axh2[i] = axh2[i-1] + fac * xj[i];
		ayh2[i] = ayh2[i-1] + fac * yj[i];
		azh2[i] = azh2[i-1] + fac * zj[i];
	}
}

static void getacch_ah3(const size_t nbod, const float mass[NPLMAX], const float xh[NPLMAX], const float yh[NPLMAX], const float zh[NPLMAX], float axh3[NPLMAX], float ayh3[NPLMAX], float azh3[NPLMAX])
{
	for(size_t i = 0; i < nbod; ++i)
	{
		axh3[i] = 0.0f;
		ayh3[i] = 0.0f;
		azh3[i] = 0.0f;
	}

	for(size_t i = 1; i < nbod - 1; ++i)
	{
		for(size_t j = i + 1; j < nbod; ++j)
		{
			float dx = xh[j] - xh[i];
			float dy = yh[j] - yh[i];
			float dz = zh[j] - zh[i];
			float rji2 = dx*dx + dy*dy + dz*dz;

			float irij3 = 1.0f / (rji2 * sqrt(rji2));
			float faci = mass[i] * irij3;
			float facj = mass[j] * irij3;

			axh3[j] -= faci*dx;
			ayh3[j] -= faci*dy;
			azh3[j] -= faci*dz;

			axh3[i] += facj*dx;
			ayh3[i] += facj*dy;
			azh3[i] += facj*dz;
		}
	}
}

void getacch_sse(const size_t nbod, __m128 mass[NPLMAX], const __m128 xj[NPLMAX], const __m128 yj[NPLMAX], const __m128 zj[NPLMAX], 
	const __m128 xh[NPLMAX], const __m128 yh[NPLMAX], const __m128 zh[NPLMAX], __m128 axh[NPLMAX], __m128 ayh[NPLMAX], __m128 azh[NPLMAX])
{
	const __m128 zeroes = _mm_set1_ps(0.0f);

	getacch_ir3_sse(nbod, xj, yj, zj, ir3j_v);
	getacch_ir3_sse(nbod, xh, yh, zh, ir3h_v);

	__m128 axh0 = zeroes, ayh0 = zeroes, azh0 = zeroes;

	getacch_ah0_sse(nbod, mass, xh, yh, zh, ir3h_v, axh0, ayh0, azh0);

	getacch_ah1_sse(nbod, mass, xh, yh, zh, xj, yj, zj, ir3h_v, ir3j_v, axh1_v, ayh1_v, azh1_v);

	getacch_ah2_sse(nbod, mass, xj, yj, zj, ir3j_v, axh2_v, ayh2_v, azh2_v);

	getacch_ah3_sse(nbod, mass, xh, yh, zh, axh3_v, ayh3_v, azh3_v);

	for(size_t i = 0; i < nbod; ++i)
	{
		axh[i] = _mm_add_ps(axh0, _mm_add_ps(axh1_v[i], _mm_add_ps(axh2_v[i], axh3_v[i])));
		ayh[i] = _mm_add_ps(ayh0, _mm_add_ps(ayh1_v[i], _mm_add_ps(ayh2_v[i], ayh3_v[i])));
		azh[i] = _mm_add_ps(azh0, _mm_add_ps(azh1_v[i], _mm_add_ps(azh2_v[i], azh3_v[i])));
	}

	//for(size_t i = 1; i < nbod; ++i)
	//{
	//	axh[i] = axh0 + axh1[i] + axh2[i] + axh3[i];
	//	ayh[i] = ayh0 + ayh1[i] + ayh2[i] + ayh3[i];
	//	azh[i] = azh0 + azh1[i] + azh2[i] + azh3[i];
	//}

}

static void getacch_ir3_sse(const size_t nbod, const __m128 x[NPLMAX], const __m128 y[NPLMAX], const __m128 z[NPLMAX], __m128 ir3[NPLMAX])
{
	const __m128 ones = _mm_set1_ps(1.0f);

	for(size_t i = 0; i < nbod; ++i)
	{
		__m128 x2 = _mm_mul_ps(x[i], x[i]);
		__m128 y2 = _mm_mul_ps(y[i], y[i]);
		__m128 z2 = _mm_mul_ps(z[i], z[i]);

		__m128 r2 = _mm_add_ps(x2, _mm_add_ps(z2, y2));

		__m128 r = _mm_sqrt_ps(r2);
		__m128 r3 = _mm_mul_ps(r2, r);
		ir3[i] = _mm_div_ps(ones, r3);
	}
}

static void getacch_ah0_sse(const size_t nbod, const __m128 mass[NPLMAX], const __m128 xh[NPLMAX], const __m128 yh[NPLMAX], const __m128 zh[NPLMAX], 
	const __m128 ir3h[NPLMAX], __m128 &axh0, __m128 &ayh0, __m128 &azh0)
{
	const __m128 zeroes = _mm_set1_ps(0.0f);
	const unsigned int MM_ALIGN16 mask[4] = { 0, 0, 0xFFFFFFFF, 0xFFFFFFFF };
	const __m128 remove_sun_and_p1 =_mm_load_ps((float*) mask);

	axh0 = zeroes;
	ayh0 = zeroes;
	azh0 = zeroes;

	__m128 ir3h0 = _mm_and_ps(ir3h[0], remove_sun_and_p1);
	__m128 fac = _mm_mul_ps(mass[0], ir3h0);

	__m128 xfac = _mm_mul_ps(fac, xh[0]);
	__m128 yfac = _mm_mul_ps(fac, yh[0]);
	__m128 zfac = _mm_mul_ps(fac, zh[0]);

	for(size_t i = 1; i < nbod; ++i)
	{
		fac = _mm_mul_ps(mass[i], ir3h[i]);
		xfac = _mm_add_ps(xfac, _mm_mul_ps(fac, xh[i]));
		yfac = _mm_add_ps(yfac, _mm_mul_ps(fac, yh[i]));
		zfac = _mm_add_ps(zfac, _mm_mul_ps(fac, zh[i]));
	}

	//float MM_ALIGN16 mx[4], mxtot;
	//_mm_store_ps(mx, xfac);
	//mxtot = mx[0] + mx[1] + mx[2] + mx[3];
	//axh0 = _mm_set1_ps(0.0f - mxtot);

	//float MM_ALIGN16 my[4], mytot;
	//_mm_store_ps(my, yfac);
	//mytot = my[0] + my[1] + my[2] + my[3];
	//ayh0 = _mm_set1_ps(0.0f - mytot);

	//float MM_ALIGN16 mz[4], mztot;
	//_mm_store_ps(mz, zfac);
	//mztot = mz[0] + mz[1] + mz[2] + mz[3];
	//azh0 = _mm_set1_ps(0.0f - mztot);

	__m128 t1a = _mm_movelh_ps(yfac, xfac);
	__m128 t1b = _mm_movehl_ps(xfac, yfac);
	__m128 t1c = _mm_movelh_ps(zeroes, zfac);
	__m128 t1d = _mm_movehl_ps(zfac, zeroes);

	__m128 t2a = _mm_add_ps(t1a, t1b);
	__m128 t2b = _mm_add_ps(t1c, t1d);

	__m128 t3a = _mm_shuffle_ps(t2a, t2b, _MM_SHUFFLE(0, 2, 0, 2));
	__m128 t3b = _mm_shuffle_ps(t2a, t2b, _MM_SHUFFLE(1, 3, 1, 3));

	__m128 tot = _mm_add_ps(t3a, t3b);

	axh0 = _mm_shuffle_ps(tot, tot, _MM_SHUFFLE(0, 0, 0, 0));
	axh0 = _mm_sub_ps(zeroes, axh0);
	ayh0 = _mm_shuffle_ps(tot, tot, _MM_SHUFFLE(1, 1, 1, 1));
	ayh0 = _mm_sub_ps(zeroes, ayh0);
	azh0 = _mm_shuffle_ps(tot, tot, _MM_SHUFFLE(2, 2, 2, 2));
	azh0 = _mm_sub_ps(zeroes, azh0);
}

static void getacch_ah1_sse(const size_t nbod, const __m128 mass[NPLMAX], const __m128 xh[NPLMAX], const __m128 yh[NPLMAX], const __m128 zh[NPLMAX], const __m128 xj[NPLMAX], const __m128 yj[NPLMAX], const __m128 zj[NPLMAX], 
	const __m128 ir3h[NPLMAX], const __m128 ir3j[NPLMAX], __m128 axh1[NPLMAX], __m128 ayh1[NPLMAX], __m128 azh1[NPLMAX])
{
	float *mass_f = (float*) mass;
	__m128 mass_sun = _mm_set1_ps(mass_f[0]);

	for(size_t i = 0; i < nbod; ++i)
	{
		__m128 ah1j = _mm_mul_ps(xj[i], ir3j[i]);
		__m128 ah1h = _mm_mul_ps(xh[i], ir3h[i]);
		__m128 diff_ah1jh = _mm_sub_ps(ah1j, ah1h);
		axh1[i] = _mm_mul_ps(mass_sun, diff_ah1jh);

		ah1j = _mm_mul_ps(yj[i], ir3j[i]);
		ah1h = _mm_mul_ps(yh[i], ir3h[i]);
		diff_ah1jh = _mm_sub_ps(ah1j, ah1h);
		ayh1[i] = _mm_mul_ps(mass_sun, diff_ah1jh);

		ah1j = _mm_mul_ps(zj[i], ir3j[i]);
		ah1h = _mm_mul_ps(zh[i], ir3h[i]);
		diff_ah1jh = _mm_sub_ps(ah1j, ah1h);
		azh1[i] = _mm_mul_ps(mass_sun, diff_ah1jh);
	}
}

static void getacch_ah2_sse(const size_t nbod, const __m128 mass[NPLMAX], const __m128 xj[NPLMAX], const __m128 yj[NPLMAX], const __m128 zj[NPLMAX], const __m128 ir3j[NPLMAX], 
	__m128 axh2[NPLMAX], __m128 ayh2[NPLMAX], __m128 azh2[NPLMAX])
{
	float *etaj_s;
	__m128 *etaj;
	__m128 mass_sun;

	etaj_s = (float*) _mm_malloc(NPLMAX * sizeof(float), 16);
	const float *mass_s = (float*)mass;

	etaj = (__m128*) etaj_s;

	etaj_s[0] = 0.0f;
	etaj_s[1] = mass_s[0];
	for(size_t i = 2; i < nbod * 4; ++i)
	{
		etaj_s[i] = etaj_s[i-1] + mass_s[i-1];
	}

	mass_sun = _mm_set1_ps(mass_s[0]);

	const unsigned int MM_ALIGN16 mask[4] = { 0, 0, 0xFFFFFFFF, 0xFFFFFFFF };
	const __m128 remove_sun_and_p1 =_mm_load_ps((float*) mask);

	for(size_t i = 0; i < nbod; ++i)
	{
		__m128 fac = _mm_div_ps(_mm_mul_ps(_mm_mul_ps(mass[i], mass_sun), ir3j[i]), etaj[i]);
		axh2[i] = _mm_mul_ps(fac, xj[i]);
		ayh2[i] = _mm_mul_ps(fac, yj[i]);
		azh2[i] = _mm_mul_ps(fac, zj[i]);
	}

	axh2[0] = _mm_and_ps(axh2[0], remove_sun_and_p1);
	ayh2[0] = _mm_and_ps(ayh2[0], remove_sun_and_p1);
	azh2[0] = _mm_and_ps(azh2[0], remove_sun_and_p1);

	//float *axh2_s = (float*)axh2;
	//float *ayh2_s = (float*)ayh2;
	//float *azh2_s = (float*)azh2;

	for(size_t i = 2; i < nbod * 4; ++i)
	{
		axh2_s[i] += axh2_s[i-1];
		ayh2_s[i] += ayh2_s[i-1];
		azh2_s[i] += azh2_s[i-1];
	}
}

static void getacch_ah3_sse(const size_t nbod, __m128 mass[NPLMAX], const __m128 xh[NPLMAX], const __m128 yh[NPLMAX], const __m128 zh[NPLMAX], __m128 axh3[NPLMAX], __m128 ayh3[NPLMAX], __m128 azh3[NPLMAX])
{
	__m128 zeroes = _mm_set1_ps(0.0f);
	__m128 ones = _mm_set1_ps(1.0f);

	for(size_t i = 0; i < nbod; ++i)
	{
		axh3[i] = zeroes;
		ayh3[i] = zeroes;
		azh3[i] = zeroes;
	}

	const unsigned int MM_ALIGN16 mask[4] = { 0, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF };
	const __m128 remove_sun =_mm_load_ps((float*) mask);

	mass[0] = _mm_and_ps(mass[0], remove_sun);

	for(size_t i = 0; i < nbod; ++i)
	{
		__m128 xacci = zeroes, yacci = zeroes, zacci = zeroes;
		__m128 xhi = xh[i], yhi = yh[i], zhi = zh[i], massi = mass[i];

		__m128 dx, dy, dz, rji2, irji3, facj, faci;

		for(size_t j = i + 1; j < nbod; ++j)
		{
			dx = _mm_sub_ps(xh[j], xhi);
			dy = _mm_sub_ps(yh[j], yhi);
			dz = _mm_sub_ps(zh[j], zhi);

			rji2 = _mm_add_ps(_mm_add_ps(_mm_mul_ps(dx, dx), _mm_mul_ps(dy, dy)), _mm_mul_ps(dz, dz));
		
			irji3 = _mm_div_ps(ones, _mm_mul_ps(rji2, _mm_sqrt_ps(rji2)));
			faci = _mm_mul_ps(massi, irji3);
			facj = _mm_mul_ps(mass[j], irji3);

			axh3[j] = _mm_sub_ps(axh3[j], _mm_mul_ps(faci, dx));
			ayh3[j] = _mm_sub_ps(ayh3[j], _mm_mul_ps(faci, dy));
			azh3[j] = _mm_sub_ps(azh3[j], _mm_mul_ps(faci, dz));

			xacci = _mm_add_ps(xacci, _mm_mul_ps(facj, dx));
			yacci = _mm_add_ps(yacci, _mm_mul_ps(facj, dy));
			zacci = _mm_add_ps(zacci, _mm_mul_ps(facj, dz));
		}

		xacci = _mm_shuffle_ps(xacci, xacci, _MM_SHUFFLE(0, 3, 2, 1));
		yacci = _mm_shuffle_ps(yacci, yacci, _MM_SHUFFLE(0, 3, 2, 1));
		zacci = _mm_shuffle_ps(zacci, zacci, _MM_SHUFFLE(0, 3, 2, 1));

		__m128 xhi_rot1 = _mm_shuffle_ps(xh[i], xh[i], _MM_SHUFFLE(0, 3, 2, 1));
		__m128 yhi_rot1 = _mm_shuffle_ps(yh[i], yh[i], _MM_SHUFFLE(0, 3, 2, 1));
		__m128 zhi_rot1 = _mm_shuffle_ps(zh[i], zh[i], _MM_SHUFFLE(0, 3, 2, 1));
		__m128 mass_rot1 = _mm_shuffle_ps(mass[i], mass[i], _MM_SHUFFLE(0, 3, 2, 1));

		dx = _mm_sub_ps(xhi_rot1, xh[i]);
		dy = _mm_sub_ps(yhi_rot1, yh[i]);
		dz = _mm_sub_ps(zhi_rot1, zh[i]);

		rji2 = _mm_add_ps(_mm_add_ps(_mm_mul_ps(dx, dx), _mm_mul_ps(dy, dy)), _mm_mul_ps(dz, dz));
		
		irji3 = _mm_div_ps(ones, _mm_mul_ps(rji2, _mm_sqrt_ps(rji2)));
		facj = _mm_mul_ps(mass_rot1, irji3);

		axh3[i] = _mm_add_ps(axh3[i], _mm_mul_ps(facj, dx));
		ayh3[i] = _mm_add_ps(ayh3[i], _mm_mul_ps(facj, dy));
		azh3[i] = _mm_add_ps(azh3[i], _mm_mul_ps(facj, dz));

		xhi = xhi_rot1;
		yhi = yhi_rot1;
		zhi = zhi_rot1;
		massi = mass_rot1;

		for(size_t j = i + 1; j < nbod; ++j)
		{
			dx = _mm_sub_ps(xh[j], xhi);
			dy = _mm_sub_ps(yh[j], yhi);
			dz = _mm_sub_ps(zh[j], zhi);

			rji2 = _mm_add_ps(_mm_add_ps(_mm_mul_ps(dx, dx), _mm_mul_ps(dy, dy)), _mm_mul_ps(dz, dz));
		
			irji3 = _mm_div_ps(ones, _mm_mul_ps(rji2, _mm_sqrt_ps(rji2)));
			faci = _mm_mul_ps(massi, irji3);
			facj = _mm_mul_ps(mass[j], irji3);

			axh3[j] = _mm_sub_ps(axh3[j], _mm_mul_ps(faci, dx));
			ayh3[j] = _mm_sub_ps(ayh3[j], _mm_mul_ps(faci, dy));
			azh3[j] = _mm_sub_ps(azh3[j], _mm_mul_ps(faci, dz));

			xacci = _mm_add_ps(xacci, _mm_mul_ps(facj, dx));
			yacci = _mm_add_ps(yacci, _mm_mul_ps(facj, dy));
			zacci = _mm_add_ps(zacci, _mm_mul_ps(facj, dz));
		}

		xacci = _mm_shuffle_ps(xacci, xacci, _MM_SHUFFLE(0, 3, 2, 1));
		yacci = _mm_shuffle_ps(yacci, yacci, _MM_SHUFFLE(0, 3, 2, 1));
		zacci = _mm_shuffle_ps(zacci, zacci, _MM_SHUFFLE(0, 3, 2, 1));

		__m128 xhi_rot2 = _mm_shuffle_ps(xh[i], xh[i], _MM_SHUFFLE(1, 0, 3, 2));
		__m128 yhi_rot2 = _mm_shuffle_ps(yh[i], yh[i], _MM_SHUFFLE(1, 0, 3, 2));
		__m128 zhi_rot2 = _mm_shuffle_ps(zh[i], zh[i], _MM_SHUFFLE(1, 0, 3, 2));
		__m128 mass_rot2 = _mm_shuffle_ps(mass[i], mass[i], _MM_SHUFFLE(1, 0, 3, 2));

		dx = _mm_sub_ps(xhi_rot2, xh[i]);
		dy = _mm_sub_ps(yhi_rot2, yh[i]);
		dz = _mm_sub_ps(zhi_rot2, zh[i]);

		rji2 = _mm_add_ps(_mm_add_ps(_mm_mul_ps(dx, dx), _mm_mul_ps(dy, dy)), _mm_mul_ps(dz, dz));
		
		irji3 = _mm_div_ps(ones, _mm_mul_ps(rji2, _mm_sqrt_ps(rji2)));
		facj = _mm_mul_ps(mass_rot2, irji3);

		axh3[i] = _mm_add_ps(axh3[i], _mm_mul_ps(facj, dx));
		ayh3[i] = _mm_add_ps(ayh3[i], _mm_mul_ps(facj, dy));
		azh3[i] = _mm_add_ps(azh3[i], _mm_mul_ps(facj, dz));

		xhi = xhi_rot2;
		yhi = yhi_rot2;
		zhi = zhi_rot2;
		massi = mass_rot2;

		for(size_t j = i + 1; j < nbod; ++j)
		{
			dx = _mm_sub_ps(xh[j], xhi);
			dy = _mm_sub_ps(yh[j], yhi);
			dz = _mm_sub_ps(zh[j], zhi);

			rji2 = _mm_add_ps(_mm_add_ps(_mm_mul_ps(dx, dx), _mm_mul_ps(dy, dy)), _mm_mul_ps(dz, dz));
		
			irji3 = _mm_div_ps(ones, _mm_mul_ps(rji2, _mm_sqrt_ps(rji2)));
			faci = _mm_mul_ps(massi, irji3);
			facj = _mm_mul_ps(mass[j], irji3);

			axh3[j] = _mm_sub_ps(axh3[j], _mm_mul_ps(faci, dx));
			ayh3[j] = _mm_sub_ps(ayh3[j], _mm_mul_ps(faci, dy));
			azh3[j] = _mm_sub_ps(azh3[j], _mm_mul_ps(faci, dz));

			xacci = _mm_add_ps(xacci, _mm_mul_ps(facj, dx));
			yacci = _mm_add_ps(yacci, _mm_mul_ps(facj, dy));
			zacci = _mm_add_ps(zacci, _mm_mul_ps(facj, dz));
		}

		xacci = _mm_shuffle_ps(xacci, xacci, _MM_SHUFFLE(0, 3, 2, 1));
		yacci = _mm_shuffle_ps(yacci, yacci, _MM_SHUFFLE(0, 3, 2, 1));
		zacci = _mm_shuffle_ps(zacci, zacci, _MM_SHUFFLE(0, 3, 2, 1));


		__m128 xhi_rot3 = _mm_shuffle_ps(xh[i], xh[i], _MM_SHUFFLE(2, 1, 0, 3));
		__m128 yhi_rot3 = _mm_shuffle_ps(yh[i], yh[i], _MM_SHUFFLE(2, 1, 0, 3));
		__m128 zhi_rot3 = _mm_shuffle_ps(zh[i], zh[i], _MM_SHUFFLE(2, 1, 0, 3));
		__m128 mass_rot3 = _mm_shuffle_ps(mass[i], mass[i], _MM_SHUFFLE(2, 1, 0, 3));

		dx = _mm_sub_ps(xhi_rot3, xh[i]);
		dy = _mm_sub_ps(yhi_rot3, yh[i]);
		dz = _mm_sub_ps(zhi_rot3, zh[i]);

		rji2 = _mm_add_ps(_mm_add_ps(_mm_mul_ps(dx, dx), _mm_mul_ps(dy, dy)), _mm_mul_ps(dz, dz));
		
		irji3 = _mm_div_ps(ones, _mm_mul_ps(rji2, _mm_sqrt_ps(rji2)));
		facj = _mm_mul_ps(mass_rot3, irji3);

		axh3[i] = _mm_add_ps(axh3[i], _mm_mul_ps(facj, dx));
		ayh3[i] = _mm_add_ps(ayh3[i], _mm_mul_ps(facj, dy));
		azh3[i] = _mm_add_ps(azh3[i], _mm_mul_ps(facj, dz));

		xhi = xhi_rot3;
		yhi = yhi_rot3;
		zhi = zhi_rot3;
		massi = mass_rot3;

		for(size_t j = i + 1; j < nbod; ++j)
		{
			dx = _mm_sub_ps(xh[j], xhi);
			dy = _mm_sub_ps(yh[j], yhi);
			dz = _mm_sub_ps(zh[j], zhi);

			rji2 = _mm_add_ps(_mm_add_ps(_mm_mul_ps(dx, dx), _mm_mul_ps(dy, dy)), _mm_mul_ps(dz, dz));
		
			irji3 = _mm_div_ps(ones, _mm_mul_ps(rji2, _mm_sqrt_ps(rji2)));
			faci = _mm_mul_ps(massi, irji3);
			facj = _mm_mul_ps(mass[j], irji3);

			axh3[j] = _mm_sub_ps(axh3[j], _mm_mul_ps(faci, dx));
			ayh3[j] = _mm_sub_ps(ayh3[j], _mm_mul_ps(faci, dy));
			azh3[j] = _mm_sub_ps(azh3[j], _mm_mul_ps(faci, dz));

			xacci = _mm_add_ps(xacci, _mm_mul_ps(facj, dx));
			yacci = _mm_add_ps(yacci, _mm_mul_ps(facj, dy));
			zacci = _mm_add_ps(zacci, _mm_mul_ps(facj, dz));
		}

		xacci = _mm_shuffle_ps(xacci, xacci, _MM_SHUFFLE(0, 3, 2, 1));
		yacci = _mm_shuffle_ps(yacci, yacci, _MM_SHUFFLE(0, 3, 2, 1));
		zacci = _mm_shuffle_ps(zacci, zacci, _MM_SHUFFLE(0, 3, 2, 1));

		axh3[i] = _mm_add_ps(axh3[i], xacci);
		ayh3[i] = _mm_add_ps(ayh3[i], yacci);
		azh3[i] = _mm_add_ps(azh3[i], zacci);

	}
}

void getacch_ah3_tp(const size_t nbod, const size_t ntp, const float mass[NPLMAX], const float xh[NPLMAX], const float yh[NPLMAX], const float zh[NPLMAX], 
	const float xht[NPLMAX], const float yht[NPLMAX], const float zht[NPLMAX], float axh3[NPLMAX], float ayh3[NPLMAX], float azh3[NPLMAX], const float mcent)
{
	for(size_t i = 0; i < ntp; ++i)
	{
		axh3[i] = 0.0f;
		ayh3[i] = 0.0f;
		azh3[i] = 0.0f;
	}


	#pragma omp parallel for
	for(size_t j = 0; j < ntp; ++j)
	{
		float rji2, irij3, fac;

		for(size_t i = 0; i < nbod; ++i)
		{
			float dx = xht[j] - xh[i];
			float dy = yht[j] - yh[i];
			float dz = zht[j] - zh[i];
			rji2 = dx*dx + dy*dy + dz*dz;

			irij3 = 1.0f / (rji2 * sqrt(rji2));
			fac = mass[i] * irij3;

			axh3[j] = axh3[j] - fac * dx;
			ayh3[j] = ayh3[j] - fac * dy;
			azh3[j] = azh3[j] - fac * dz;
		}

		rji2 = xht[j]*xht[j] + yht[j]*yht[j] + zht[j]*zht[j];
		irij3 = 1.0f / (rji2 * sqrt(rji2));

		fac = mcent*irij3;
		axh3[j] = axh3[j] + fac * xht[j];
		ayh3[j] = ayh3[j] + fac * yht[j];
		azh3[j] = azh3[j] + fac * zht[j];
	}
}

void getacch_tp(const size_t nbod, const size_t npl, const size_t ntp, const float mass[NPLMAX], const float xh[NPLMAX], const float yh[NPLMAX], const float zh[NPLMAX], 
	const float xht[NPLMAX], const float yht[NPLMAX], const float zht[NPLMAX], float axht[NPLMAX], float ayht[NPLMAX], float azht[NPLMAX], const bool bar, const float mcent)
{
	//float ir3ht[NPLMAX];
	float ir3h[NPLMAX];

	getacch_ir3(nbod, 1, xh, yh, zh, ir3h);
	//getacch_ir3(ntp, 0, xht, yht, zht, ir3ht);

	float axh0, ayh0, azh0;

	getacch_ah0(1, nbod, mass, xh, yh, zh, ir3h, axh0, ayh0, azh0);

	//printf("axh0, ayh0, azh0: %E, %E, %E\n", axh0, ayh0, azh0);

	float axh3[NTPMAX];
	float ayh3[NTPMAX];
	float azh3[NTPMAX];

	//printf("mcent: %E\n", mcent);

	getacch_ah3_tp(nbod, ntp, mass, xh, yh, zh, xht, yht, zht, axh3, ayh3, azh3, mcent);

	//#pragma omp parallel for
	for(size_t i = 0; i < ntp; ++i)
	{
		axht[i] = axh0 + axh3[i];
		//printf("%E, %E, %E\n", axh3[i], axh0, axht[i]);
		ayht[i] = ayh0 + ayh3[i];
		azht[i] = azh0 + azh3[i];
	}
}