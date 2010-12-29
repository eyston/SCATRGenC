#include "specifics.h"
#include <cstdio> // use C++ headers, not C (when they exist).
#include <omp.h>
#include <xmmintrin.h>

#include "structures.h"
#include "io_input.h"
#include "getacch.h"
#include "getacch_sse.h"
#include "coord.h"

int main(int argc, char* argv[]) {

	// can't really use std::auto_ptr or std::uniq_ptr because of _mm_free; simply...
	typedef storage_t<NPLMAX, NTPMAX> store_t;
	store_t * const store = allocate<store_t>();
	store_t::scalar_t &scalar = store->scalar; // save some typing.
	store_t::sse_t &sse = store->sse;


	size_t ntp = xio_input_particles("tpverybig.in", scalar.particles);

	// too damn annoying.
	float * __restrict const mass_s = scalar.planets.MASS.m;
	#define PLA3(X) scalar.planets.X.x, scalar.planets.X.y, scalar.planets.X.z
	#define PAR3(X) scalar.particles.X.x, scalar.particles.X.y, scalar.particles.X.z
	// #define PAR1(X) scalar.planets.X.m

	{
		size_t nbod, npl;
		xio_input_planets("pl.in.8", nbod, npl, scalar.planets);
		xcoord_h2j(nbod, scalar.planets);
		double start = omp_get_wtime();
		for(size_t i = 0; i < 1000; ++i)
			getacch_tp(nbod, npl, ntp, mass_s, PLA3(H), PAR3(HT), PAR3(AHT), false, mass_s[0]);


		double end = omp_get_wtime();
		printf("getacch_tp time: %f\n", end - start);
	}


	{
		size_t nbod, npl;
		xio_input_planets("pl.in.8", nbod, npl, scalar.planets);
		xcoord_h2j(nbod, scalar.planets);

		double start = omp_get_wtime();
		for(size_t i = 0; i < 1000000; ++i)
			// why 2?
			xgetacch_sse(2, sse.planets);
		double end = omp_get_wtime();
		printf("getacch_sse time: %f\n", end - start);
	}

	{
		size_t nbod, npl;
		xio_input_planets("pl.in.8", nbod, npl, scalar.planets);
		xcoord_h2j(nbod, scalar.planets);

		double start = omp_get_wtime();
		for(size_t i = 0; i < 1000000; ++i)
			getacch(nbod, mass_s, PLA3(J), PLA3(H), PLA3(AH));

		double end = omp_get_wtime();
		printf("getacch time: %f\n", end - start);
	}


	// you were leaking.
	deallocate(store);
	return 0;

}