#include <stdio.h>
#include "structures.h"

void coord_h2j(const size_t nbod, const float mass[NPLMAX], const float xh[NPLMAX], const float yh[NPLMAX], const float zh[NPLMAX], const float vxh[NPLMAX], const float vyh[NPLMAX], const float vzh[NPLMAX],
	float xj[NPLMAX], float yj[NPLMAX], float zj[NPLMAX], float vxj[NPLMAX], float vyj[NPLMAX], float vzj[NPLMAX])
{
	float eta[NPLMAX];
	eta[0] = mass[0];

	for(size_t i = 1; i < nbod; ++i)
	{
		eta[i] = eta[i-1] + mass[i];
	}

	xj[0] = 0.0f;
	yj[0] = 0.0f;
	zj[0] = 0.0f;
	vxj[0] = 0.0f;
	vyj[0] = 0.0f;
	vzj[0] = 0.0f;

	xj[1] = xh[1];
	yj[1] = yh[1];
	zj[1] = zh[1];
	vxj[1] = vxh[1];
	vyj[1] = vyh[1];
	vzj[1] = vzh[1];

	float sumx = mass[1] * xh[1];
	float sumy = mass[1] * yh[1];
	float sumz = mass[1] * zh[1];
	float sumvx = mass[1] * vxh[1];
	float sumvy = mass[1] * vyh[1];
	float sumvz = mass[1] * vzh[1];

	float capx = sumx / eta[1];
	float capy = sumy / eta[1];
	float capz = sumz / eta[1];
	float capvx = sumvx / eta[1];
	float capvy = sumvy / eta[1];
	float capvz = sumvz / eta[1];

	for(size_t i = 2; i < nbod; ++i)
	{
		xj[i] = xh[i] - capx;
		yj[i] = yh[i] - capy;
		zj[i] = zh[i] - capz;
		vxj[i] = vxh[i] - capx;
		vyj[i] = vyh[i] - capy;
		vzj[i] = vzh[i] - capz;

		if(i < nbod - 1)
		{
			sumx += mass[i] * xh[i];
			sumy += mass[i] * yh[i];
			sumz += mass[i] * zh[i];
			sumvx += mass[i] * vxh[i];
			sumvy += mass[i] * vyh[i];
			sumvz += mass[i] * vzh[i];

			capx = sumx/eta[i];
			capy = sumy/eta[i];
			capz = sumz/eta[i];
			capvx = sumvx/eta[i];
			capvy = sumvy/eta[i];
			capvz = sumvz/eta[i];
		}
	}
}

void xcoord_h2j(const size_t nbod, planets_scalar_t &planets) {
	// looks like we're precomputing something.
	// doesn't matter. not gonna re-write.
	const float *mass = &planets.MASS[0];
	const float *xh = planets.H.x, *yh = planets.H.y, *zh = planets.H.z;
	const float *vxh = planets.VH.x, *vyh = planets.VH.y, *vzh = planets.VH.z;
	float *xj = planets.J.x, *yj = planets.J.y, *zj = planets.J.z;
	float *vxj = planets.VJ.x, *vyj = planets.VJ.y, *vzj = planets.VJ.z;

	float eta[NPLMAX];
	eta[0] = mass[0];

	for(size_t i = 1; i < nbod; ++i)
	{
		eta[i] = eta[i-1] + mass[i];
	}

	xj[0] = 0.0f;
	yj[0] = 0.0f;
	zj[0] = 0.0f;
	vxj[0] = 0.0f;
	vyj[0] = 0.0f;
	vzj[0] = 0.0f;

	xj[1] = xh[1];
	yj[1] = yh[1];
	zj[1] = zh[1];
	vxj[1] = vxh[1];
	vyj[1] = vyh[1];
	vzj[1] = vzh[1];

	float sumx = mass[1] * xh[1];
	float sumy = mass[1] * yh[1];
	float sumz = mass[1] * zh[1];
	float sumvx = mass[1] * vxh[1];
	float sumvy = mass[1] * vyh[1];
	float sumvz = mass[1] * vzh[1];

	float capx = sumx / eta[1];
	float capy = sumy / eta[1];
	float capz = sumz / eta[1];
	float capvx = sumvx / eta[1];
	float capvy = sumvy / eta[1];
	float capvz = sumvz / eta[1];

	for(size_t i = 2; i < nbod; ++i)
	{
		xj[i] = xh[i] - capx;
		yj[i] = yh[i] - capy;
		zj[i] = zh[i] - capz;
		vxj[i] = vxh[i] - capx;
		vyj[i] = vyh[i] - capy;
		vzj[i] = vzh[i] - capz;

		if(i < nbod - 1)
		{
			sumx += mass[i] * xh[i];
			sumy += mass[i] * yh[i];
			sumz += mass[i] * zh[i];
			sumvx += mass[i] * vxh[i];
			sumvy += mass[i] * vyh[i];
			sumvz += mass[i] * vzh[i];

			capx = sumx/eta[i];
			capy = sumy/eta[i];
			capz = sumz/eta[i];
			capvx = sumvx/eta[i];
			capvy = sumvy/eta[i];
			capvz = sumvz/eta[i];
		}
	}
}
