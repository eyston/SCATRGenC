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


struct Vector
{
	float x;
	float y;
	float z;
};

struct NBody
{
	float mass;
	float rpl;
	
	Vector h_loc;
	Vector h_vel;
	
	Vector j_loc;
	Vector j_vel;
};

struct NBodyParameters
{
	int nbod, npl;
	float j2rp2, j4rp4;
	NBody *bodies;
};

