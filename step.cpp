#include <stdio.h>
#include "step.h"
#include "coord.h"
#include "kickvh.h"
#include "getacch.h"
#include "drift.h"

void step_kdk(bool i1st, float time, int nbod, int npl, int ntp, float *mass, float j2rp2, float j4rp4, float *xh, float *yh, float *zh, float *vxh, float *vyh, float *vzh,
	float *xht, float *yht, float *zht, float *vxht, float *vyht, float *vzht, int istat, float rstat, float dt, int bar)
{
	printf("step_kdk called\n");

	step_kdk_pl(i1st, nbod, mass, j2rp2, j4rp4, xh, yh, zh, vxh, vyh, vzh, dt);

	float *xbeg = new float[nbod];
	float *ybeg = new float[nbod];
	float *zbeg = new float[nbod];

	step_kdk_tp(i1st, time, nbod, npl, ntp, mass, j2rp2, j4rp4, xbeg, ybeg, zbeg, xh, yh, zh, xht, yht, zht, vxht, vyht, vzht, istat, dt, bar);
}

void step_kdk_pl(bool i1st, int nbod, float *mass, float j2rp2, float j4rp4, float *xh, float *yh, float *zh, float *vxh, float *vyh, float *vzh, float dt)
{
	printf("step_kdk_pl called\n");

	// coord_h2j();
	// getacch();
	// kickvh();
	// coord_vh2vj();
	// drift();
	// coord_j2h();
	// getacch();
	// kickvh();
}

void step_kdk_tp(bool i1st, float time, int nbod, int npl, int ntp, float *mass, float j2rp2, float j4rp4, float *xbeg, float *ybeg, float *zbeg, float *xend, float *yend, float *zend,
	float *xht, float *yht, float *zht, float *vxht, float *vyht, float *vzht, int istat, float dt, int bar)
{
	printf("step_kdk_tp called\n");

	// getacch_tp();
	// kickvh_tp();
	// drift_tp();
	// getacch_tp();
	// kickvh_tp();
}