/*
 * io_input.cpp
 *
 *  Created on: Nov 23, 2010
 *      Author: Huey
 */

#include <stdio.h>
#include "io_input.h"
#include "structures.h"

SimulationParameters io_input_params(const char *input_file_name)
{
	SimulationParameters setup = SimulationParameters();

	FILE *input_file;
	input_file = fopen(input_file_name, "r");

	fscanf(input_file, "%f %f %f %f\n", &setup.t0, &setup.tstop, &setup.dt, &setup.tinc);
	fscanf(input_file, "%f %f\n", &setup.dtout, &setup.dtdump);
	fscanf(input_file, "%c %c %c %c %c %c\n", &setup.lflag_c[0], &setup.lflag_c[1], &setup.lflag_c[2], &setup.lflag_c[3], &setup.lflag_c[4], &setup.lflag_c[5]);
	fscanf(input_file, "%f %f %f %f %c\n", &setup.rmin, &setup.rmax, &setup.rmaxu, &setup.qmin, &setup.lclose_c);
	fscanf(input_file, "%f\n", &setup.rcrit);
	fscanf(input_file, "%s\n", setup.output_file_name);
	fscanf(input_file, "%s\n", setup.frame);
	fscanf(input_file, "%s\n", setup.fopenstat);

	fclose(input_file);

	printf("%f, %f, %f, %f\n", setup.t0, setup.tstop, setup.dt, setup.tinc);
	printf("%f, %f\n", setup.dtout, setup.dtdump);
	printf("%c, %c, %c, %c, %c, %c\n", setup.lflag_c[0], setup.lflag_c[1], setup.lflag_c[2], setup.lflag_c[3], setup.lflag_c[4], setup.lflag_c[5]);

	for(int i=0; i<6; i++)
	{
		setup.lflag[i] = setup.lflag_c[i] == 'T';
	}

	printf("%d, %d, %d, %d, %d, %d\n", setup.lflag[0], setup.lflag[1], setup.lflag[2], setup.lflag[3], setup.lflag[4], setup.lflag[5]);

	setup.lclose = setup.lclose_c == 'T';

	printf("%f, %f, %f, %f, %d\n", setup.rmin, setup.rmax, setup.rmaxu, setup.qmin, setup.lclose);
	printf("%f\n", setup.rcrit);
	printf("%s\n", setup.output_file_name);
	printf("%s\n", setup.frame);
	printf("%s\n", setup.fopenstat);

	return setup;
}

NBodies io_input_planets(const char *input_file_name)
{
	FILE *input_file;
	input_file = fopen(input_file_name, "r");

	int nbod, npl;

	fscanf(input_file, "%i %i\n", &nbod, &npl);
	printf("%i, %i\n", nbod, npl);

	NBodies nbodies = NBodies();
	nbodies.SetNBod(nbod);
	nbodies.npl = npl;

	fscanf(input_file, "%f %f %f\n", &nbodies.mass[0], &nbodies.j2rp2, &nbodies.j4rp4);
	fscanf(input_file, "%f %f %f\n", &nbodies.xh[0], &nbodies.yh[0], &nbodies.zh[0]);
	fscanf(input_file, "%f %f %f\n", &nbodies.vxh[0], &nbodies.vyh[0], &nbodies.vzh[0]);

	for(int i=1; i < nbodies.nbod; i++)
	{
		fscanf(input_file, "%f %f\n", &nbodies.mass[i], &nbodies.rpl[i]);
		fscanf(input_file, "%f %f %f\n", &nbodies.xh[i], &nbodies.yh[i], &nbodies.zh[i]);
		fscanf(input_file, "%f %f %f\n", &nbodies.vxh[i], &nbodies.vyh[i], &nbodies.vzh[i]);
	}

	fclose(input_file);

	return nbodies;
}
