#include <stdio.h>
#include "structures.h"
#include "rmvs3.h"
#include "rmvs3inner.h"

void rmvs3_step(bool i1st, float t, int nbod, int npl, int ntp, float *mass, float j2rp2, float j4rp4, float *xh, float *yh, float *zh, float *vxh, float *vyh, float *vzh,
	float *xht, float *yht, float *zht, float *vxht, float *vyht, float *vzht, int istat, int rstat, float dt, float tinc, float fopenstat, int nclose, int iclose,
	float rmin, float rmax, float rmaxu, float qmin, float rplsq, float rcrit)
{
	printf("rmvs3_step called\n");

	int lclose = 0, partmap = 0;

	rmvs3_inner(i1st, t, nbod, npl, ntp, mass, j2rp2, j4rp4, xh, yh, zh, vxh, vyh, vzh, xht, yht, zht, vxht, vyht, vzht, istat, rstat, dt, tinc, fopenstat, partmap,
		lclose, rmin, rmax, rmaxu, qmin, rplsq);

	// coord_h2b();

	// coord_h2b();

	// rmvs3_chk();

	// step_kdk_tp();

	// rmvs3_interp();

	// rmvs3_step_out();

	// step_kdk_tp();
}

void rmvs3_chk(int npl, int nbod, int ntp, float *mass, float *xh, float *yh, float *zh, float *vxh, float *vyh, float *vzh,
	float *xht, float *yht, float *zht, float *vxht, float *vyht, float *vzht, int istat, float dt, float rts, int icflg, int nenc, int itpenc, int ienc, int istart, int passhill)
{
	printf("rmvs3_chk called\n");

	//util_hills();

	//rmvs_chk_ind();
}

void rmvs3_interp(int nbod, float *xbeg, float *ybeg, float *zbeg, float *vxbeg, float *vybeg, float *vzbeg, float *xend, float *yend, float *zend, float *vxend, float *vyend, float *vzend,
	float dt, float msun, int nt, float *xtmp, float *ytmp, float *ztmp, float *vxtmp, float *vytmp, float *vztmp, int istart)
{
	printf("rmvs3_interp\n");

	// not sure what this does
}