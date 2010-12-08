#include <stdio.h>
#include <omp.h>

#include "structures.h"
#include "io_input.h"
#include "getacch.h"
#include "coord.h"


int main(int argc, char* argv[]) {

	SimulationParameters params = io_input_params("param.in");

	NBodies bodies = io_input_planets("pl.in");

	size_t nbod = bodies.nbod;
	float j2rp2 = bodies.j2rp2;
	float j4rp4 = bodies.j4rp4;
	float xj[NPLMAX], yj[NPLMAX], zj[NPLMAX], xh[NPLMAX], yh[NPLMAX], zh[NPLMAX];
	float vxj[NPLMAX], vyj[NPLMAX], vzj[NPLMAX], vxh[NPLMAX], vyh[NPLMAX], vzh[NPLMAX];
	float mass[NPLMAX], axh[NPLMAX], ayh[NPLMAX], azh[NPLMAX];

	for(size_t i = 0; i < nbod; ++i)
	{
		mass[i] = bodies.mass[i];
		xh[i] = bodies.xh[i];
		yh[i] = bodies.yh[i];
		zh[i] = bodies.zh[i];
		vxh[i] = bodies.vxh[i];
		vyh[i] = bodies.vyh[i];
		vzh[i] = bodies.vzh[i];
	}

	coord_h2j(nbod, mass, xh, yh, zh, vxh, vyh, vzh, xj, yj, zj, vxj, vyj, vzj);

	double start = omp_get_wtime();

	for(size_t i = 0; i < 10000000; ++i)
	{
		getacch(nbod, mass, j2rp2, j4rp4, xj, yj, zj, xh, yh, zh, axh, ayh, azh);
	}

	double end = omp_get_wtime();

	//getacch(nbod, mass, j2rp2, j4rp4, xj, yj, zj, xh, yh, zh, axh, ayh, azh);

	printf("%E\n", mass[1]);
	printf("%E, %E, %E\n", xh[1], yh[1], zh[1]);
	printf("%E, %E, %E\n", vxh[1], vyh[1], vzh[1]);
	printf("%E, %E, %E\n", xj[1], yj[1], zj[1]);
	printf("%E, %E, %E\n", vxj[1], vyj[1], vzj[1]);
	printf("%E, %E, %E\n", axh[1], ayh[1], azh[1]);

	printf("time: %f\n", end - start);

	return 0;
}
