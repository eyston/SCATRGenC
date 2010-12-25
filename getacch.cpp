#include <stdio.h>
#include <xmmintrin.h>
#include <math.h>
#include "structures.h"
#include "getacch.h"

void getacch(const size_t nbod, const float* __restrict mass, const float* __restrict xj, const float* __restrict yj, const float* __restrict zj, 
	const float* __restrict xh, const float* __restrict yh, const float* __restrict zh, float* __restrict axh, float* __restrict ayh, float* __restrict azh)
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

static void getacch_ir3(const size_t nbod, const size_t istart, const float* __restrict x, const float* __restrict y, const float* __restrict z, float* __restrict ir3)
{
	for(size_t i = istart; i < nbod; ++i)
	{
		float r2 = x[i]*x[i] + y[i]*y[i] + z[i]*z[i];
		float ir = 1.0f / sqrt(r2);
		ir3[i] = ir / r2;
	}
}

static void getacch_ah0(const size_t istart, const size_t nbod, const float* __restrict mass, const float* __restrict xh, const float* __restrict yh, const float* __restrict zh, 
	const float* __restrict ir3h, float &axh0, float &ayh0, float &azh0)
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

static void getacch_ah1(const size_t nbod, const float* __restrict mass, const float* __restrict xh, const float* __restrict yh, const float* __restrict zh, const float* __restrict xj, const float* __restrict yj, const float* __restrict zj, 
	const float* __restrict ir3h, const float* __restrict ir3j, float* __restrict axh1, float* __restrict ayh1, float* __restrict azh1)
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

static void getacch_ah2(const size_t nbod, const float* __restrict mass, const float* __restrict xj, const float* __restrict yj, const float* __restrict zj, const float* __restrict ir3j, 
	float* __restrict axh2, float* __restrict ayh2, float* __restrict azh2)
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

static void getacch_ah3(const size_t nbod, const float* __restrict mass, const float* __restrict xh, const float* __restrict yh, const float* __restrict zh, float* __restrict axh3, float* __restrict ayh3, float* __restrict azh3)
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

void getacch_ah3_tp(const size_t nbod, const size_t ntp, const float* __restrict mass, const float* __restrict xh, const float* __restrict yh, const float* __restrict zh, 
	const float* __restrict xht, const float* __restrict yht, const float* __restrict zht, float* __restrict axh3, float* __restrict ayh3, float* __restrict azh3, const float mcent)
{
	for(size_t i = 0; i < ntp; ++i)
	{
		axh3[i] = 0.0f;
		ayh3[i] = 0.0f;
		azh3[i] = 0.0f;
	}


	//#pragma omp parallel for
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

void getacch_tp(const size_t nbod, const size_t npl, const size_t ntp, const float* __restrict mass, const float* __restrict xh, const float* __restrict yh, const float* __restrict zh, 
	const float* __restrict xht, const float* __restrict yht, const float* __restrict zht, float* __restrict axht, float* __restrict ayht, float* __restrict azht, const bool bar, const float mcent)
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