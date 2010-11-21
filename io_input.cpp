struct SimulationSetup
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

SimulationSetup io_input(const char *input_file_name)
{
	SimulationSetup setup = SimulationSetup();

	FILE *input_file;
	input_file = fopen("param.in", "r");

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
