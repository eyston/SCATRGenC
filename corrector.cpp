#include <math.h>
#include "coord_h2j.h"

void corectterm_hel(int nbod, NBody bodies[])
{
	coord_h2j(nbod, bodies);
}

void invcorrector(int nbod, NBody bodies[], float tinc, float dt)
{
	
	/* float *xh2 = new float[nbod];
	float *yh2 = new float[nbod];
	float *zh2 = new float[nbod];
	float *vxh2 = new float[nbod];
	float *vyh2 = new float[nbod];
	float *vzh2 = new float[nbod];
	
	float gamma, a1, a2, b1, b2, adt, bdt;
	
	
	for(int i=0; i < nbod; i++)
	{
		xh2[i] = bodies[i].xh;
		yh2[i] = bodies[i].yh;
		zh2[i] = bodies[i].zh;
		vxh2[i] = bodies[i].vxh;
		vyh2[i] = bodies[i].vyh;
		vzh2[i] = bodies[i].vzh;
	}
	
	gamma = sqrt(10.0);
	a2 = gamma / 0.5;
	b2 = gamma / 24.0;
	adt = dt * a2 / tinc;
	bdt = dt * b2 / tinc;
	
	// call correctterm_hel
	
	a1 = 3.0 * gamma / 10.0;
	b1 = 0.0 - gamma / 72.0;
	adt = dt * a1 / tinc;
	bdt = dt * b1 / tinc;
	
	// call correctterm_hel
	
	adt = dt * a2;
	bdt = dt * b2;
	
	// call correctterm_bar
	
	adt = dt * a1;
	bdt = dt * b1;

	// call correctterm_bar */

}