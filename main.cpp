#include <iostream>
#include <fstream>
#include <stdio.h>

#include "structures.h"
#include "coord_h2j.h"
#include "io_input.cpp"
#include "corrector.cpp"

using namespace std;


int main()
{

	SimulationParameters params = io_input_params("param.in");
	
	NBodyParameters pl_params = io_input_planets("pl.in");
	
	int t = params.t0;
	int tout = params.t0 + params.dtout;
	int tdump = params.t0 + params.dtdump;
	
	invcorrector(pl_params.nbod, pl_params.bodies, params.tinc, params.dt);
	
	return 0;
}