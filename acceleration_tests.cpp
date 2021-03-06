#include <limits.h>
#include <gtest/gtest.h>
#include <cmath>
#include <vector>

#include "io_input.h"
#include "coord.h"
#include "structures.h"
#include "specifics.h"
#include "getacch_sse.h"
#include "getacch.h"

#define ASSERT_NEAR_PERCENT(EXPECTED, ACTUAL, PERCENT) ASSERT_NEAR(EXPECTED, ACTUAL, fabs(EXPECTED * PERCENT));

typedef storage_t<NPLMAX, NTPMAX> store_t;

class Inner_Accelerations_Fortran : public virtual testing::Test
{
	protected:
		Inner_Accelerations_Fortran(const char * input_file_name) : sun_offset(1)
		{
			file_name = input_file_name;
		}

		~Inner_Accelerations_Fortran()
		{
		}

		void read_in_fortran_ir3j()
		{
			parse_values_for_lines_that_start_with("ir3j", fortran_ir3j);
		}

		void read_in_fortran_ir3h()
		{
			parse_values_for_lines_that_start_with("ir3h", fortran_ir3h);
		}

		void read_in_fortran_ah0()
		{
			vec3_scalar_t ah0;

			parse_values_for_lines_that_start_with("axh0", ah0);
			fortran_axh0 = ah0.x[sun_offset];
			fortran_ayh0 = ah0.y[sun_offset];
			fortran_azh0 = ah0.z[sun_offset];
		}

		void read_in_fortran_ah1()
		{
			parse_values_for_lines_that_start_with("axh1", fortran_ah1);
		}

		void read_in_fortran_ah2()
		{
			parse_values_for_lines_that_start_with("axh2", fortran_ah2);
		}

		void read_in_fortran_ah3()
		{
			parse_values_for_lines_that_start_with("axh3", fortran_ah3);
		}

	vec1_scalar_t fortran_ir3j;
	vec1_scalar_t fortran_ir3h;
	vec3_scalar_t fortran_ah1;
	vec3_scalar_t fortran_ah2;
	vec3_scalar_t fortran_ah3;
	float fortran_axh0, fortran_ayh0, fortran_azh0;

	private:

		enum { length = 1024 };

		void parse_values_for_lines_that_start_with(const char* identifier, vec1_scalar_t &values)
		{
			parse_values_for_lines_that_start_with(identifier, "%*s %*s %f", values);
		}

		void parse_values_for_lines_that_start_with(const char* identifier, const char * format_string, vec1_scalar_t &values)
		{
			std::vector<std::string> lines = all_lines_that_start_with(identifier);

			for(size_t i = 0; i < lines.size(); ++i)
			{
				float value;
				sscanf(lines[i].c_str(), format_string, &value);
				values.m[i + sun_offset] = value;
			}
		}

		void parse_values_for_lines_that_start_with(const char* identifier, vec3_scalar_t &values)
		{
			parse_values_for_lines_that_start_with(identifier, "%*s %*s %*s %f %f %f", values);
		}

		void parse_values_for_lines_that_start_with(const char* identifier, const char * format_string, vec3_scalar_t &values)
		{
			std::vector<std::string> lines = all_lines_that_start_with(identifier);

			for(size_t i = 0; i < lines.size(); ++i)
			{
				float value_x, value_y, value_z;
				sscanf(lines[i].c_str(), format_string, &value_x, &value_y, &value_z);
				values.x[i + sun_offset] = value_x;
				values.y[i + sun_offset] = value_y;
				values.z[i + sun_offset] = value_z;
			}
		}

		std::vector<std::string> all_lines_that_start_with(const char* identifier)
		{
			std::vector<std::string> lines;
			std::string identifier_with_comma = std::string(identifier).append(",");

			if (FILE *input_file = fopen(file_name, "r")) {
				char line[length];
				char string[length] = "";

				while(fgets(line, length, input_file) != NULL)
				{
					sscanf(line, "%s", string);
					if(identifier_with_comma.compare(string) == 0)
					{
						lines.push_back(line);
					}
				}

				fclose(input_file);
			}

			return lines;
		}


	const char * file_name;
	const int sun_offset;
};

class Inner_Accelerations_Fortran_8_Planets : public Inner_Accelerations_Fortran
{
	protected:
		Inner_Accelerations_Fortran_8_Planets() : Inner_Accelerations_Fortran("results.8.out") { }
};

// sanity check tests to make sure we are reading files correctly
// ultimately I want results.X.out where X is 4, 5, 6, 7, 8, 9

TEST_F(Inner_Accelerations_Fortran_8_Planets, should_match_file_values_for_single_floats)
{
	read_in_fortran_ah0();

	ASSERT_NEAR_PERCENT(-4.29461297470980564E-010f, fortran_axh0, 0.00001f);
	ASSERT_NEAR_PERCENT(1.17541738894071035E-010f, fortran_ayh0, 0.00001f);
	ASSERT_NEAR_PERCENT(-6.93754178176094684E-010f, fortran_azh0, 0.00001f);
}

TEST_F(Inner_Accelerations_Fortran_8_Planets, should_match_file_values_for_vector_1)
{
	read_in_fortran_ir3j();

	ASSERT_NEAR_PERCENT(8.10678918941581603E-003f, fortran_ir3j.m[1], 0.00001f);
	ASSERT_NEAR_PERCENT(9.89038081836128395E-004f, fortran_ir3j.m[2], 0.00001f);
	ASSERT_NEAR_PERCENT(1.63144591141573295E-004f, fortran_ir3j.m[3], 0.00001f);
	ASSERT_NEAR_PERCENT(3.66552844843908994E-005f, fortran_ir3j.m[4], 0.00001f);
	ASSERT_NEAR_PERCENT(4.76525394977701514E-019f, fortran_ir3j.m[5], 0.00001f);
	ASSERT_NEAR_PERCENT(4.17659140665433338E-018f, fortran_ir3j.m[6], 0.00001f);
	ASSERT_NEAR_PERCENT(7.92853122035865894E-018f, fortran_ir3j.m[7], 0.00001f);
}

TEST_F(Inner_Accelerations_Fortran_8_Planets, should_match_file_values_for_vector_3)
{
	read_in_fortran_ah1();

	ASSERT_NEAR_PERCENT(0.0000000000000000f, fortran_ah1.x[1], 0.00001f);
	ASSERT_NEAR_PERCENT(0.0000000000000000f, fortran_ah1.y[1], 0.00001f);
	ASSERT_NEAR_PERCENT(0.0000000000000000f, fortran_ah1.z[1], 0.00001f);

	ASSERT_NEAR_PERCENT(1.86521048010972244E-010f, fortran_ah1.x[2], 0.00001f);
	ASSERT_NEAR_PERCENT(-1.41709727243463246E-009f, fortran_ah1.y[2], 0.00001f);
	ASSERT_NEAR_PERCENT(3.09476539847044786E-010f, fortran_ah1.z[2], 0.00001f);

	ASSERT_NEAR_PERCENT(5.34878335975960556E-011f, fortran_ah1.x[3], 0.00001f);
	ASSERT_NEAR_PERCENT(5.34305085351076559E-010f, fortran_ah1.y[3], 0.00001f);
	ASSERT_NEAR_PERCENT(1.12973779135292565E-010f, fortran_ah1.z[3], 0.00001f);

	ASSERT_NEAR_PERCENT(-2.60399717363032296E-011, fortran_ah1.x[4], 0.00001f);
	ASSERT_NEAR_PERCENT(6.62848716577778652E-011f, fortran_ah1.y[4], 0.00001f);
	ASSERT_NEAR_PERCENT(-4.21779742715585840E-011f, fortran_ah1.z[4], 0.00001f);

	ASSERT_NEAR_PERCENT(6.97796363932148213E-025f, fortran_ah1.x[5], 0.00001f);
	ASSERT_NEAR_PERCENT(1.38418267816305365E-024f, fortran_ah1.y[5], 0.00001f);
	ASSERT_NEAR_PERCENT(-4.17774939020519457E-025f, fortran_ah1.z[5], 0.00001f);

	ASSERT_NEAR_PERCENT(-7.37423237304377330E-017f, fortran_ah1.x[6], 0.00001f);
	ASSERT_NEAR_PERCENT(-5.54409854623788595E-017f, fortran_ah1.y[6], 0.00001f);
	ASSERT_NEAR_PERCENT(7.67904399881583287E-017f, fortran_ah1.z[6], 0.00001f);

	ASSERT_NEAR_PERCENT(2.25585231396107693E-016f, fortran_ah1.x[7], 0.00001f);
	ASSERT_NEAR_PERCENT(2.94631809056399085E-016f, fortran_ah1.y[7], 0.00001f);
	ASSERT_NEAR_PERCENT(3.49540464733356542E-017f, fortran_ah1.z[7], 0.00001f);
}

union MM_ALIGN64 accelerations_t {
	struct scalar_t {
		vec1_scalar_t ir3j;
		vec1_scalar_t ir3h;
		vec3_scalar_t ah1;
		vec3_scalar_t ah2;
		vec3_scalar_t ah3;
	};
	struct sse_t {
		vec1_sse_t ir3j;
		vec1_sse_t ir3h;
		vec3_sse_t ah1;
		vec3_sse_t ah2;
		vec3_sse_t ah3;
	};

	scalar_t scalar;
	sse_t sse;
};

class Inner_Accelerations_Vector : public virtual testing::Test
{
	protected:
		Inner_Accelerations_Vector(const char * input_file_name, const int nbod) :
				accelerations(allocate<accelerations_t>()),
				vector_ir3j(accelerations->scalar.ir3j),
				vector_ir3h(accelerations->scalar.ir3h),
				vector_ah1(accelerations->scalar.ah1),
				vector_ah2(accelerations->scalar.ah2),
				vector_ah3(accelerations->scalar.ah3),
				store(allocate<store_t>()),
				nbod_v(nbod)
		{
			file_name = input_file_name;
		}

		~Inner_Accelerations_Vector()
		{
		}

		void calculate_vector_ir3j()
		{
			initialize_planets();
			x_getacch_ir3_sse_test(nbod_v, store->sse.planets.J, accelerations->sse.ir3j);
		}

		void calculate_vector_ir3h()
		{
			initialize_planets();
			x_getacch_ir3_sse_test(nbod_v, store->sse.planets.H, accelerations->sse.ir3h);
		}

		void calculate_vector_ah0()
		{
			initialize_planets();
			calculate_vector_ir3h();
			x_getacch_ah0_sse_test(nbod_v, store->sse.planets.MASS, accelerations->sse.ir3h, store->sse.planets.H, vector_axh0, vector_ayh0, vector_azh0);
		}

		void calculate_vector_ah1()
		{
			initialize_planets();
			calculate_vector_ir3h();
			calculate_vector_ir3j();
			x_getacch_ah1_sse_test(nbod_v, store->sse.planets.MASS, accelerations->sse.ir3h, accelerations->sse.ir3j, store->sse.planets.H, store->sse.planets.J, accelerations->sse.ah1);
		}

		void calculate_vector_ah2()
		{
			initialize_planets();
			calculate_vector_ir3j();
			x_getacch_ah2_sse_test(nbod_v, store->sse.planets.MASS, accelerations->sse.ir3j, store->sse.planets.J, accelerations->sse.ah2);
		}

		void calculate_vector_ah3()
		{
			initialize_planets();
			x_getacch_ah3_sse_test(nbod_v, store->sse.planets.MASS, store->sse.planets.H, accelerations->sse.ah3);
		}

	private:

		void initialize_planets()
		{
			size_t nbod, npl;
			xio_input_planets(file_name, nbod, npl, store->scalar.planets);
			xcoord_h2j(nbod, store->scalar.planets);
		}

	accelerations_t *accelerations;

	protected:

	vec1_scalar_t &vector_ir3j;
	vec1_scalar_t &vector_ir3h;
	float vector_axh0, vector_ayh0, vector_azh0;
	vec3_scalar_t &vector_ah1;
	vec3_scalar_t &vector_ah2;
	vec3_scalar_t &vector_ah3;

	private:

	store_t * store;

	const char * file_name;
	const int nbod_v;
};

class Inner_Accelerations_Vector_8_Planets : public Inner_Accelerations_Vector
{
	protected:
		Inner_Accelerations_Vector_8_Planets() : Inner_Accelerations_Vector("pl.in.8", 2) { }
};

TEST_F(Inner_Accelerations_Vector_8_Planets, should_match_file_values_for_single_floats)
{
	calculate_vector_ir3j();

	ASSERT_NEAR_PERCENT(8.10678918941581603E-003f, vector_ir3j[1], 0.00001f);
}

class Inner_Accelerations_Vector_And_Fortran_8_Planets : public Inner_Accelerations_Vector_8_Planets, public Inner_Accelerations_Fortran_8_Planets
{
protected:
	Inner_Accelerations_Vector_And_Fortran_8_Planets() : Inner_Accelerations_Vector_8_Planets(), Inner_Accelerations_Fortran_8_Planets() { }
};

TEST_F(Inner_Accelerations_Vector_And_Fortran_8_Planets, should_match_fortran_and_vector_values_for_ir3j)
{
	read_in_fortran_ir3j();
	calculate_vector_ir3j();

	for(int i = 1; i < 8; ++i)
	{
		ASSERT_NEAR_PERCENT(fortran_ir3j[i], vector_ir3j[i], 0.00001f);
	}
}

TEST_F(Inner_Accelerations_Vector_And_Fortran_8_Planets, should_match_fortran_and_vector_values_for_ir3h)
{
	read_in_fortran_ir3h();
	calculate_vector_ir3h();

	for(int i = 1; i < 8; ++i)
	{
		ASSERT_NEAR_PERCENT(fortran_ir3h[i], vector_ir3h[i], 0.00001f);
	}
}

TEST_F(Inner_Accelerations_Vector_And_Fortran_8_Planets, should_match_fortran_and_vector_values_for_ah0)
{
	read_in_fortran_ah0();
	calculate_vector_ah0();

	ASSERT_NEAR_PERCENT(fortran_axh0, vector_axh0, 0.00001f);
	ASSERT_NEAR_PERCENT(fortran_ayh0, vector_ayh0, 0.00001f);
	ASSERT_NEAR_PERCENT(fortran_azh0, vector_azh0, 0.00001f);
}


// E: tolerances are bad on this one ... we need to verify code!
//TEST_F(Inner_Accelerations_Vector_And_Fortran_8_Planets, should_match_fortran_and_vector_values_for_ah1)
//{
//	read_in_fortran_ah1();
//	calculate_vector_ah1();
//
//	for(int i = 1; i < 8; ++i)
//	{
//		ASSERT_NEAR_PERCENT(fortran_ah1.x[i], vector_ah1.x[i], 0.00001f);
//		ASSERT_NEAR_PERCENT(fortran_ah1.y[i], vector_ah1.y[i], 0.00001f);
//		ASSERT_NEAR_PERCENT(fortran_ah1.z[i], vector_ah1.z[i], 0.00001f);
//	}
//}

TEST_F(Inner_Accelerations_Vector_And_Fortran_8_Planets, should_match_fortran_and_vector_values_for_ah2)
{
	read_in_fortran_ah2();
	calculate_vector_ah2();

	for(int i = 1; i < 8; ++i)
	{
		ASSERT_NEAR_PERCENT(fortran_ah2.x[i], vector_ah2.x[i], 0.00001f);
		ASSERT_NEAR_PERCENT(fortran_ah2.y[i], vector_ah2.y[i], 0.00001f);
		ASSERT_NEAR_PERCENT(fortran_ah2.z[i], vector_ah2.z[i], 0.00001f);
	}
}

TEST_F(Inner_Accelerations_Vector_And_Fortran_8_Planets, should_match_fortran_and_vector_values_for_ah3)
{
	read_in_fortran_ah3();
	calculate_vector_ah3();

	for(int i = 1; i < 8; ++i)
	{
		ASSERT_NEAR_PERCENT(fortran_ah3.x[i], vector_ah3.x[i], 0.00001f);
		ASSERT_NEAR_PERCENT(fortran_ah3.y[i], vector_ah3.y[i], 0.00001f);
		ASSERT_NEAR_PERCENT(fortran_ah3.z[i], vector_ah3.z[i], 0.00001f);
	}
}
