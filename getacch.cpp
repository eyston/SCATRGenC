#include <stdio.h>
#include <math.h>
#include "structures.h"
#include "getacch.h"

void getacch(const size_t nbod, const float mass[NPLMAX], const float j2rp2, const float j4rp4, const float xj[NPLMAX], const float yj[NPLMAX], const float zj[NPLMAX], 
	const float xh[NPLMAX], const float yh[NPLMAX], const float zh[NPLMAX], float axh[NPLMAX], float ayh[NPLMAX], float azh[NPLMAX])
{
	float ir3j[NPLMAX]; float irj[NPLMAX];
	float ir3h[NPLMAX]; float irh[NPLMAX];

	getacch_ir3(nbod, 1, xj, yj, zj, ir3j, irj);
	getacch_ir3(nbod, 1, xh, yh, zh, ir3h, irh);

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

void getacch_ir3(const size_t nbod, const size_t istart, const float x[NPLMAX], const float y[NPLMAX], const float z[NPLMAX], float ir3[NPLMAX], float ir[NPLMAX])
{
	for(size_t i = istart; i < nbod; ++i)
	{
		float r2 = x[i]*x[i] + y[i]*y[i] + z[i]*z[i];
		ir[i] = 1.0f / sqrt(r2);
		ir3[i] = ir[i] / r2;
	}
}

void getacch_ah0(const size_t istart, const size_t nbod, const float mass[NPLMAX], const float xh[NPLMAX], const float yh[NPLMAX], const float zh[NPLMAX], 
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

void getacch_ah1(const size_t nbod, const float mass[NPLMAX], const float xh[NPLMAX], const float yh[NPLMAX], const float zh[NPLMAX], const float xj[NPLMAX], const float yj[NPLMAX], const float zj[NPLMAX], 
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

void getacch_ah2(const size_t nbod, const float mass[NPLMAX], const float xj[NPLMAX], const float yj[NPLMAX], const float zj[NPLMAX], const float ir3j[NPLMAX], 
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
		float fac = mass[i] * mass[1] * ir3j[i] / etaj;
		axh2[i] = axh2[i-1] + fac * xj[i];
		ayh2[i] = ayh2[i-1] + fac * yj[i];
		azh2[i] = azh2[i-1] + fac * zj[i];
	}
}

void getacch_ah3(const size_t nbod, const float mass[NPLMAX], const float xh[NPLMAX], const float yh[NPLMAX], const float zh[NPLMAX], float axh3[NPLMAX], float ayh3[NPLMAX], float azh3[NPLMAX])
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