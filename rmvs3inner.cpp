#include <stdio.h>
#include "rmvs3inner.h"
#include "drift.h"
#include "step.h"

void rmvs3_hel(int nbod, int npl, int ntp, float t, float *mass, float j2rp2, float j4rp4, float *xh, float *yh, float *zh, float *vxh, float *vyh, float *vzh,
	float *xht, float *yht, float *zht, float *vxht, float *vyht, float *vzht, int istat, float dt, float tinc, int iclose, int nclose, int rcrit)
{
	printf("rmvs3_hel called\n");

	float mtot = 0.0f, xhtdum = 0.0f, yhtdum = 0.0f, zhtdum = 0.0f, vxhtdum = 0.0f, vyhtdum = 0.0f, vzhtdum = 0.0f, r2crit = 0.0f;

	int iflag = 0;

	rmvs3_chk_hel(nbod, mtot, xhtdum, yhtdum, zhtdum, vxhtdum, vyhtdum, vzhtdum, dt, iflag, r2crit);

	// corrector();

	// invcorrector();
}

void rmvs3_chk_hel(int nbod, float mtot, float xhtdum, float yhtdum, float zhtdum, float vxhtdum, float vyhtdum, float vzhtdum, float dt, int iflag, float r2crit)
{
	printf("rmvs3_chk_hel called\n");

	drift_one(mtot, xhtdum, yhtdum, zhtdum, vxhtdum, vyhtdum, vzhtdum, dt);
}

void rmvs3_inner(bool i1st, float time, int nbod, int npl, int ntp, float *mass, float j2rp2, float j4rp4, float *xh, float *yh, float *zh, float *vxh, float *vyh, float *vzh,
	float *xht, float *yht, float *zht, float *vxht, float *vyht, float *vzht, int istat, int rstat, float dt, float tinc, float fopenstat, float partmap, float lclose,
	float rmin, float rmax, float rmaxu, float qmin, float rplsq)
{
	printf("rmvs3_inner called\n");

	// rmvs3_chk();
	
	int bar = 0;
	step_kdk(i1st, time, nbod, npl, ntp, mass, j2rp2, j4rp4, xh, yh, zh, vxh, vyh, vzh, xht, yht, zht, vxht, vyht, vzht, istat, rstat, dt, bar);

	step_kdk_pl(i1st, nbod, mass, j2rp2, j4rp4, xh, yh, zh, vxh, vyh, vzh, dt);

	// rmvs3_interp();

	// rmvs3_step_out();

	float *xbeg = new float[nbod];
	float *ybeg = new float[nbod];
	float *zbeg = new float[nbod];

	step_kdk_tp(i1st, time, nbod, npl, ntp, mass, j2rp2, j4rp4, xbeg, ybeg, zbeg, xh, yh, zh, xht, yht, zht, vxht, vyht, vzht, istat, dt, bar);

	// discard();
}