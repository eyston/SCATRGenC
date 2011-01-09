#include <limits.h>
#include <gtest/gtest.h>
#include <cmath>

#include "io_input.h"
#include "coord.h"
#include "structures.h"
#include "specifics.h"
#include "getacch_sse.h"
#include "getacch.h"

#define ASSERT_NEAR_PERCENT(EXPECTED, ACTUAL, PERCENT) ASSERT_NEAR(EXPECTED, ACTUAL, fabs(EXPECTED * PERCENT));

typedef storage_t<NPLMAX, NTPMAX> store_t;


class Acceleration_8_FortranPlanetSimulation : public testing::Test
{
	protected:
		Acceleration_8_FortranPlanetSimulation() : 
			store_v(allocate<store_t>()), 
			sse(store_v->sse),
			vector_accelerations(store_v->scalar.planets.AH) { }

		~Acceleration_8_FortranPlanetSimulation()
		{
			deallocate(store_v);
		}

		void readInFortranAccelerations()
		{
			enum { length = 256 };
			if (FILE *input_file = fopen("results.8.out", "r")) {
				char line[length];
				char string[length] = "";
				char string2[length] = "";
				char string3[length] = "";

				while(strcmp(string, "axh,") != 0)
				{
					fgets(line, length, input_file);
					sscanf(line, "%s", string);
				}

				int i = 1;

				while(strcmp(string, "axh,") == 0)
				{
					float axh, ayh, azh;
					sscanf(line, "%s %s %s %f %f %f", string, string2, string3, &axh, &ayh, &azh);
					fortran_accelerations.x[i] = axh;
					fortran_accelerations.y[i] = ayh;
					fortran_accelerations.z[i] = azh;

					i++;

					if(fgets(line, length, input_file) == NULL)
						break;
				}

				fclose(input_file);
			}
		}

		void calculateVectorAccelerations()
		{
			#define PLA3V(X) sse.planets.X.x, sse.planets.X.y, sse.planets.X.z

			{
				size_t nbod, npl;
				xio_input_planets("pl.in.8", nbod, npl, store_v->scalar.planets);
				xcoord_h2j(nbod, store_v->scalar.planets);

				xgetacch_sse(2, sse.planets);
				//xgetacch_sse(2, sse.planets);
				//getacch_sse(2, sse.planets.MASS.m, PLA3V(J), PLA3V(H), PLA3V(AH));
			}
		}

	store_t * store_v;

	store_t::sse_t &sse;

	vec3_scalar_t &vector_accelerations;
	vec3_scalar_t fortran_accelerations;
};

TEST_F(Acceleration_8_FortranPlanetSimulation, VectorAndScalarPlanetAccelerationsMatch)
{
	readInFortranAccelerations();
	calculateVectorAccelerations();

	for(size_t i = 1; i < 8; ++i)
	{
		ASSERT_NEAR_PERCENT(fortran_accelerations.x[i], vector_accelerations.x[i], 0.001f);
		ASSERT_NEAR_PERCENT(fortran_accelerations.y[i], vector_accelerations.y[i], 0.001f);
		ASSERT_NEAR_PERCENT(fortran_accelerations.z[i], vector_accelerations.z[i], 0.001f);
	}
}

class Acceleration_8_PlanetSimulation : public testing::Test
{
	protected:
		Acceleration_8_PlanetSimulation() : 
			store_s(allocate<store_t>()), 
			store_v(allocate<store_t>()), 
			scalar(store_s->scalar),
			sse(store_v->sse),
			scalar_accelerations(store_s->scalar.planets.AH),
			vector_accelerations(store_v->scalar.planets.AH) { }

		~Acceleration_8_PlanetSimulation()
		{
			deallocate(store_s);
			deallocate(store_v);
		}


		void calculateScalarAccelerations()
		{
			float * __restrict const mass_s = scalar.planets.MASS.m;
			#define PLA3(X) scalar.planets.X.x, scalar.planets.X.y, scalar.planets.X.z
			#define PAR3(X) scalar.particles.X.x, scalar.particles.X.y, scalar.particles.X.z

			{
				size_t nbod, npl;
				xio_input_planets("pl.in.8", nbod, npl, scalar.planets);
				xcoord_h2j(nbod, scalar.planets);
				getacch(nbod, mass_s, PLA3(J), PLA3(H), PLA3(AH));
			}
		}

		void calculateVectorAccelerations()
		{
			#define PLA3V(X) sse.planets.X.x, sse.planets.X.y, sse.planets.X.z

			{
				size_t nbod, npl;
				xio_input_planets("pl.in.8", nbod, npl, store_v->scalar.planets);
				xcoord_h2j(nbod, store_v->scalar.planets);

				xgetacch_sse(2, sse.planets);
				//xgetacch_sse(2, sse.planets);
				//getacch_sse(2, sse.planets.MASS.m, PLA3V(J), PLA3V(H), PLA3V(AH));
			}
		}

	store_t * store_s;
	store_t * store_v;

	store_t::scalar_t &scalar;
	store_t::sse_t &sse;

	vec3_scalar_t &scalar_accelerations;
	vec3_scalar_t &vector_accelerations;
};


TEST_F(Acceleration_8_PlanetSimulation, VectorAndScalarPlanetAccelerationsMatch)
{
	calculateScalarAccelerations();
	calculateVectorAccelerations();

	for(size_t i = 1; i < scalar_accelerations.capacity; ++i)
	{
		ASSERT_NEAR_PERCENT(scalar_accelerations.x[i], vector_accelerations.x[i], 0.001f);
		ASSERT_NEAR_PERCENT(scalar_accelerations.y[i], vector_accelerations.y[i], 0.001f);
		ASSERT_NEAR_PERCENT(scalar_accelerations.z[i], vector_accelerations.z[i], 0.001f);
	}
}

//class Acceleration_9_PlanetSimulation : public testing::Test
//{
//	protected:
//		Acceleration_9_PlanetSimulation() : 
//			store_s(allocate<store_t>()), 
//			store_v(allocate<store_t>()), 
//			scalar(store_s->scalar),
//			sse(store_v->sse),
//			scalar_accelerations(store_s->scalar.planets.AH),
//			vector_accelerations(store_v->scalar.planets.AH) { }
//
//		~Acceleration_9_PlanetSimulation()
//		{
//			deallocate(store_s);
//			deallocate(store_v);
//		}
//
//
//		void calculateScalarAccelerations()
//		{
//			float * __restrict const mass_s = scalar.planets.MASS.m;
//			#define PLA3(X) scalar.planets.X.x, scalar.planets.X.y, scalar.planets.X.z
//			#define PAR3(X) scalar.particles.X.x, scalar.particles.X.y, scalar.particles.X.z
//
//			{
//				size_t nbod, npl;
//				xio_input_planets("pl.in.9", nbod, npl, scalar.planets);
//				xcoord_h2j(nbod, scalar.planets);
//				getacch(nbod, mass_s, PLA3(J), PLA3(H), PLA3(AH));
//			}
//		}
//
//		void calculateVectorAccelerations()
//		{
//			#define PLA3V(X) sse.planets.X.x, sse.planets.X.y, sse.planets.X.z
//
//			{
//				size_t nbod, npl;
//
//				sun_t<__m128> sun;
//
//				yio_input_planets("pl.in.9", nbod, npl, sun, store_v->scalar.planets);
//				//xcoord_h2j(nbod, store_v->scalar.planets);
//
//				//ygetacch_sse(2, sun, sse.planets);
//			}
//		}
//
//	store_t * store_s;
//	store_t * store_v;
//
//	store_t::scalar_t &scalar;
//	store_t::sse_t &sse;
//
//	vec3_scalar_t &scalar_accelerations;
//	vec3_scalar_t &vector_accelerations;
//};


//TEST_F(Acceleration_9_PlanetSimulation, VectorAndScalarPlanetAccelerationsMatch)
//{
//	calculateScalarAccelerations();
//	calculateVectorAccelerations();
//
//	for(size_t i = 1; i < 8; ++i)
//	{
//		ASSERT_NEAR_PERCENT(scalar_accelerations.x[i], vector_accelerations.x[i], 0.001f);
//		ASSERT_NEAR_PERCENT(scalar_accelerations.y[i], vector_accelerations.y[i], 0.001f);
//		ASSERT_NEAR_PERCENT(scalar_accelerations.z[i], vector_accelerations.z[i], 0.001f);
//	}
//}