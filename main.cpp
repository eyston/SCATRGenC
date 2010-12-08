#include <stdio.h>

#include "structures.h"
#include "io_input.h"
#include "corrector.h"

int main(int argc, char* argv[]) {

	SimulationParameters params = io_input_params("param.in");

	NBodies bodies = io_input_planets("pl.in");

//	int t = params.t0;
//	int tout = params.t0 + params.dtout;
//	int tdump = params.t0 + params.dtdump;

	invcorrector(bodies, params.tinc, params.dt);

	return 0;
}
