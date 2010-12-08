#ifndef GETACHH_H_
#define GETACHH_H_

void getacch(const size_t nbod, const float mass[NPLMAX], const float j2rp2, const float j4rp4, const float xj[NPLMAX], const float yj[NPLMAX], const float zj[NPLMAX], 
	const float xh[NPLMAX], const float yh[NPLMAX], const float zh[NPLMAX], float axh[NPLMAX], float ayh[NPLMAX], float azh[NPLMAX]);

void getacch_ir3(const size_t nbod, const size_t istart, const float x[NPLMAX], const float y[NPLMAX], const float z[NPLMAX], float ir3[NPLMAX], float ir[NPLMAX]);

void getacch_ah0(const size_t istart, const size_t nbod, const float mass[NPLMAX], const float xh[NPLMAX], const float yh[NPLMAX], const float zh[NPLMAX], 
	const float ir3h[NPLMAX], float &axh0, float &ayh0, float &azh0);

void getacch_ah1(const size_t nbod, const float mass[NPLMAX], const float xh[NPLMAX], const float yh[NPLMAX], const float zh[NPLMAX], const float xj[NPLMAX], const float yj[NPLMAX], const float zj[NPLMAX], 
	const float ir3h[NPLMAX], const float ir3j[NPLMAX], float axh1[NPLMAX], float ayh1[NPLMAX], float azh1[NPLMAX]);

void getacch_ah2(const size_t nbod, const float mass[NPLMAX], const float xj[NPLMAX], const float yj[NPLMAX], const float zj[NPLMAX], const float ir3j[NPLMAX], 
	float axh2[NPLMAX], float ayh2[NPLMAX], float azh2[NPLMAX]);

void getacch_ah3(const size_t nbod, const float mass[NPLMAX], const float xh[NPLMAX], const float yh[NPLMAX], const float zh[NPLMAX], float axh3[NPLMAX], float ayh3[NPLMAX], float azh3[NPLMAX]);

#endif /* GETACHH_H_ */
