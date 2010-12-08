/*
 * getachh.cpp
 *
 *  Created on: Nov 23, 2010
 *      Author: Huey
 */

#include <stdio.h>
#include "structures.h"
#include "getacch.h"

void getacch(NBodies bodies, float j2rp2, float j4rp4)
{
	printf("getachh called\n");

	float *ir3j = new float[bodies.nbod];
	float *irj = new float[bodies.nbod];

	getacch_ir3(bodies.nbod, 2, bodies.xj, bodies.yj, bodies.zj, ir3j, irj);

	float *ir3h = new float[bodies.nbod];
	float *irh = new float[bodies.nbod];

	getacch_ir3(bodies.nbod, 2, bodies.xh, bodies.yh, bodies.zh, ir3h, irh);

	//getacch_ah0();

	//getacch_ah1();

	//getacch_ah2();

	//getacch_ah3();

	//obl_acc();
}

void getacch_ir3(int nbod, int istart, float *x, float *y, float *z, float *ir3, float *ir)
{
	printf("getacch_ir3 called\n");
}
