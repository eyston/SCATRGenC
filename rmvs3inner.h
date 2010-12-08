#ifndef RMVS3INNER_H_
#define RMVS3INNER_H_

void rmvs3_hel(int nbod, int npl, int ntp, float t, float *mass, float j2rp2, float j4rp4, float *xh, float *yh, float *zh, float *vxh, float *vyh, float *vzh,
	float *xht, float *yht, float *zht, float *vxht, float *vyht, float *vzht, int istat, float dt, float tinc, int iclose, int nclose, int rcrit);

void rmvs3_chk_hel(int nbod, float mtot, float xhtdum, float yhtdum, float zhtdum, float vxhtdum, float vyhtdum, float vzhtdum, float dt, int iflag, float r2crit);

void rmvs3_inner(bool i1st, float time, int nbod, int npl, int ntp, float *mass, float j2rp2, float j4rp4, float *xh, float *yh, float *zh, float *vxh, float *vyh, float *vzh,
	float *xht, float *yht, float *zht, float *vxht, float *vyht, float *vzht, int istat, int rstat, float dt, float tinc, float fopenstat, float partmap, float lclose,
	float rmin, float rmax, float rmaxu, float qmin, float rplsq);

#endif /* RMVS3INNER_H_ */