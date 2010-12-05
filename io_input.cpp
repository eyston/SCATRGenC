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

NBodyParameters io_input_planets(const char *input_file_name)
{
	FILE *input_file;
	input_file = fopen(input_file_name, "r");

	NBodyParameters pl_params = NBodyParameters();
	

	fscanf(input_file, "%i %i\n", &pl_params.nbod, &pl_params.npl);
	printf("%i, %i\n", pl_params.nbod, pl_params.npl);
	
	NBody *bodies = new NBody[pl_params.nbod];
	pl_params.bodies = bodies;
	
	fscanf(input_file, "%f %f %f\n", &bodies[0].mass, &pl_params.j2rp2, &pl_params.j4rp4);
	fscanf(input_file, "%f %f %f\n", &bodies[0].h_loc.x, &bodies[0].h_loc.y, &bodies[0].h_loc.z);
	fscanf(input_file, "%f %f %f\n", &bodies[0].h_vel.x, &bodies[0].h_vel.y, &bodies[0].h_vel.z);

	printf("%f, %f, %f\n", bodies[0].mass, pl_params.j2rp2, pl_params.j4rp4);
	printf("%f, %f, %f\n", bodies[0].h_loc.x, bodies[0].h_loc.y, bodies[0].h_loc.z);
	printf("%f, %f, %f\n", bodies[0].h_vel.x, bodies[0].h_vel.y, bodies[0].h_vel.z);
	
	for(int i=1; i < pl_params.nbod; i++)
	{
		fscanf(input_file, "%f %f\n", &bodies[i].mass, &bodies[i].rpl);
		fscanf(input_file, "%f %f %f\n", &bodies[i].h_loc.x, &bodies[i].h_loc.y, &bodies[i].h_loc.z);
		fscanf(input_file, "%f %f %f\n", &bodies[i].h_vel.x, &bodies[i].h_vel.y, &bodies[i].h_vel.z);

		printf("%f, %f\n", bodies[i].mass, bodies[i].rpl);
		printf("%f, %f, %f\n", bodies[i].h_loc.x, bodies[i].h_loc.y, bodies[i].h_loc.z);
		printf("%f, %f, %f\n", bodies[i].h_vel.x, bodies[i].h_vel.y, bodies[i].h_vel.z);		
	}
	
	fclose(input_file);
	
	return pl_params;
}
