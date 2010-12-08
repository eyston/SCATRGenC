#ifndef STRUCTURES_H_
#define STRUCTURES_H_

#include <cstdlib>

const size_t NPLMAX = 1024;
const size_t NTPMAX = 1024;

struct SimulationParameters
{
	float t0, tstop, dt, tinc;
	float dtout, dtdump;
	bool lflag[6];
	char lflag_c[6];
	float rmin, rmax, rmaxu, qmin;
	bool lclose;
	char lclose_c;
	float rcrit;
	char output_file_name[256];
	char frame[4];
	char fopenstat[8];
};


struct NBodies
{
	int nbod, npl;
	float j2rp2, j4rp4;

	float *mass;
	float *rpl;

	float *xh;
	float *yh;
	float *zh;

	float *vxh;
	float *vyh;
	float *vzh;

	float *xj;
	float *yj;
	float *zj;

	float *vxj;
	float *vyj;
	float *vzj;

	void SetNBod(int p_nbod)
	{
		nbod = p_nbod;

		mass = new float[nbod];
		rpl = new float[nbod];

		xh = new float[nbod];
		yh = new float[nbod];
		zh = new float[nbod];

		vxh = new float[nbod];
		vyh = new float[nbod];
		vzh = new float[nbod];

		xj = new float[nbod];
		yj = new float[nbod];
		zj = new float[nbod];

		vxj = new float[nbod];
		vyj = new float[nbod];
		vzj = new float[nbod];
	}
};


#endif /* STRUCTURES_H_ */
