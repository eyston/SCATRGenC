#ifndef STEP_H_
#define STEP_H_

void step_kdk(bool i1st, float time, int nbod, int npl, int ntp, float *mass, float j2rp2, float j4rp4, float *xh, float *yh, float *zh, float *vxh, float *vyh, float *vzh,
	float *xht, float *yht, float *zht, float *vxht, float *vyht, float *vzht, int istat, float rstat, float dt, int bar);

void step_kdk_pl(bool i1st, int nbod, float *mass, float j2rp2, float j4rp4, float *xh, float *yh, float *zh, float *vxh, float *vyh, float *vzh, float dt);

void step_kdk_tp(bool i1st, float time, int nbod, int npl, int ntp, float *mass, float j2rp2, float j4rp4, float *xbeg, float *ybeg, float *zbeg, float *xend, float *yend, float *zend,
	float *xht, float *yht, float *zht, float *vxht, float *vyht, float *vzht, int istat, float dt, int bar);

#endif /* STEP_H_ */