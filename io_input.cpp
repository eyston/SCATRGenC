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

void io_input_planets(const char *input_file_name, size_t &nbod, size_t &npl, float mass[NPLMAX], float rpl[NPLMAX], float xh[NPLMAX], float yh[NPLMAX], float zh[NPLMAX], float vxh[NPLMAX], float vyh[NPLMAX], float vzh[NPLMAX])
{
	const size_t length = 256;

	FILE *input_file;
	char line[length];
	input_file = fopen(input_file_name, "r");

	fgets(line, length, input_file);
	sscanf(line, "%i %i\n", &nbod, &npl);

	fgets(line, length, input_file);
	sscanf(line, "%f\n", &mass[0]);

	fgets(line, length, input_file);
	sscanf(line, "%f %f %f\n", &xh[0], &yh[0], &zh[0]);

	fgets(line, length, input_file);
	sscanf(line, "%f %f %f\n", &vxh[0], &vyh[0], &vzh[0]);

	for(size_t i=1; i < nbod; i++)
	{
		fscanf(input_file, "%f %f\n", &mass[i], &rpl[i]);
		fscanf(input_file, "%f %f %f\n", &xh[i], &yh[i], &zh[i]);
		fscanf(input_file, "%f %f %f\n", &vxh[i], &vyh[i], &vzh[i]);
	}

	fclose(input_file);
}

void io_input_particles(const char *input_file_name, size_t &ntp, float xht[NTPMAX], float yht[NTPMAX], float zht[NTPMAX], float vxht[NTPMAX], float vyht[NTPMAX], float vzht[NTPMAX])
{
	const size_t length = 256;

	FILE *input_file;
	char line[length];
	input_file = fopen(input_file_name, "r");

	fgets(line, length, input_file);
	sscanf(line, "%i", &ntp);

	for(size_t i = 0; i < ntp; i++)
	{
		fgets(line, length, input_file);
		sscanf(line, "%f %f %f", &xht[i], &yht[i], &zht[i]);
		fgets(line, length, input_file);
		sscanf(line, "%f %f %f", &vxht[i], &vyht[i], &vzht[i]);

		fgets(line, length, input_file);
		fgets(line, length, input_file);
		fgets(line, length, input_file);
		fgets(line, length, input_file);
	}

	fclose(input_file);
}

size_t xio_input_particles(const char * const input_file_name, particles_scalar_t &particles) {
	enum { length = 256 };
	size_t ntp = 0;
	if (FILE *input_file = fopen(input_file_name, "r")) {
		char line[length];
		fgets(line, length, input_file);
		sscanf(line, "%i", &ntp);

		for(size_t i = 0; i < ntp; i++) {
			fgets(line, length, input_file);
			sscanf(line, "%f %f %f", &particles.HT.x[i], &particles.HT.y[i], &particles.HT.z[i]);
			fgets(line, length, input_file);
			sscanf(line, "%f %f %f", &particles.VHT.x[i], &particles.VHT.y[i], &particles.VHT.z[i]);

			fgets(line, length, input_file);
			fgets(line, length, input_file);
			fgets(line, length, input_file);
			fgets(line, length, input_file);
		}

		fclose(input_file);
	}
	return ntp;
}

void xio_input_planets(const char * const input_file_name, size_t &nbod, size_t &npl, planets_scalar_t &planets) {
	enum { length = 256 };
	nbod = 0;
	npl = 0;
	if (FILE *input_file = fopen(input_file_name, "r")) {
		char line[length];
		fgets(line, length, input_file);
		sscanf(line, "%i %i\n", &nbod, &npl);

		// WTF?
		fgets(line, length, input_file);
		sscanf(line, "%f\n", &planets.MASS[0]);
		//FIXME: RPL[0] is never initialized or what?

		fgets(line, length, input_file);
		sscanf(line, "%f %f %f\n", &planets.H.x[0], &planets.H.y[0], &planets.H.z[0]);

		fgets(line, length, input_file);
		sscanf(line, "%f %f %f\n", &planets.VH.x[0], &planets.VH.y[0], &planets.VH.z[0]);

		for(size_t i=1; i < nbod; i++)
		{
			fscanf(input_file, "%f %f\n", &planets.MASS[i], &planets.RPL[i]);
			fscanf(input_file, "%f %f %f\n", &planets.H.x[i], &planets.H.y[i], &planets.H.z[i]);
			fscanf(input_file, "%f %f %f\n", &planets.VH.x[i], &planets.VH.y[i], &planets.VH.z[i]);
		}
		fclose(input_file);
	}
}
