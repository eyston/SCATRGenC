#ifndef RMVS3_H_
#define RMVS3_H_


void rmvs3_step(bool i1st, float t, int nbod, int npl, int ntp, float *mass, float j2rp2, float j4rp4, float *xh, float *yh, float *zh, float *vxh, float *vyh, float *vzh,
	float *xht, float *yht, float *zht, float *vxht, float *vyht, float *vzht, int istat, int rstat, float dt, float tinc, float fopenstat, int nclose, int iclose,
	float rmin, float rmax, float rmaxu, float qmin, float rplsq, float rcrit);

void rmvs3_step_out(bool i1st, int nbod, int npl, int ntp, float *mass, float j2rp2, float j4rp4, float *xbeg, float *ybeg, float *zbeg, float *vxbeg, float *vybeg, float *vzbeg,
	float *xtmp, float *ytmp, float *ztmp, float *vxtmp, float *vytmp, float *vztmp, float *xht, float *yht, float *zht, float *vxht, float *vyht, float *vzht,
	int istat, int ienc, float dt, int isperi, float peri, float time, int istart, float bar);

void rmvs3_chk(int npl, int nbod, int ntp, float *mass, float *xh, float *yh, float *zh, float *vxh, float *vyh, float *vzh,
	float *xht, float *yht, float *zht, float *vxht, float *vyht, float *vzht, int istat, float dt, float rts, int icflg, int nenc, int itpenc, int ienc, int istart, int passhill);

void rmvs3_interp(int nbod, float *xbeg, float *ybeg, float *zbeg, float *vxbeg, float *vybeg, float *vzbeg, float *xend, float *yend, float *zend, float *vxend, float *vyend, float *vzend,
	float dt, float msun, int nt, float *xtmp, float *ytmp, float *ztmp, float *vxtmp, float *vytmp, float *vztmp, int istart);

#endif /* RMVS3_H_ */