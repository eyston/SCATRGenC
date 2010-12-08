/*
 * corrector.cpp
 *
 *  Created on: Nov 23, 2010
 *      Author: Huey
 */

#include <stdio.h>
#include "structures.h"
#include "corrector.h"
#include "coord.h"
#include "drift.h"
#include "getacch.h"
#include "kickvh.h"

void invcorrector(NBodies bodies, float tinc, float dt)
{
	printf("invcorrector called\n");

	// call correctterm_hel
	correctterm_hel(bodies);

	// call correctterm_hel
	correctterm_hel(bodies);

	// call correctterm_bar
	correctterm_bar(bodies);

	// call correctterm_bar
	correctterm_bar(bodies);
}

void correctterm_hel(NBodies bodies)
{
	printf("correctterm_hel called\n");

	// call coord_h2j

	float adt;
	float j2rp2, j4rp4;

	coord_h2j(bodies);

	drift(bodies, adt);

	coord_j2h(bodies);

	getacch(bodies, j2rp2, j4rp4);

	kickvh(bodies);

	coord_h2j(bodies);

	drift(bodies, adt);

	coord_j2h(bodies);

	coord_h2j(bodies);

	drift(bodies, adt);

	coord_j2h(bodies);

	getacch(bodies, j2rp2, j4rp4);

	kickvh(bodies);

	coord_h2j(bodies);

	drift(bodies, adt);

	coord_j2h(bodies);
}

void correctterm_bar(NBodies bodies)
{
	printf("correctterm_bar called\n");

	float adt;

	coord_h2j(bodies);

	drift(bodies, adt);

	coord_j2h(bodies);

	// convert pl's to barycentric coords

	coord_h2j(bodies);

	drift(bodies, adt);

	coord_j2h(bodies);

	coord_h2j(bodies);

	drift(bodies, adt);

	coord_j2h(bodies);

	// convert pl's to barycentric coords


}
