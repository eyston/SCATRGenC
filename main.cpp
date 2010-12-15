#include <stdio.h>
#include <omp.h>
#include <xmmintrin.h>

#include "structures.h"
#include "io_input.h"
#include "getacch.h"
#include "coord.h"


int main(int argc, char* argv[]) {

	size_t nbod, npl;

	float *mass = (float*) _mm_malloc(NPLMAX * sizeof(float), 16);
	float *xh = (float*) _mm_malloc(NPLMAX * sizeof(float), 16);
	float *yh = (float*) _mm_malloc(NPLMAX * sizeof(float), 16);
	float *zh = (float*) _mm_malloc(NPLMAX * sizeof(float), 16);
	float *xj = (float*) _mm_malloc(NPLMAX * sizeof(float), 16);
	float *yj = (float*) _mm_malloc(NPLMAX * sizeof(float), 16);
	float *zj = (float*) _mm_malloc(NPLMAX * sizeof(float), 16);

	float *axh = (float*) _mm_malloc(NPLMAX * sizeof(float), 16);
	float *ayh = (float*) _mm_malloc(NPLMAX * sizeof(float), 16);
	float *azh = (float*) _mm_malloc(NPLMAX * sizeof(float), 16);

	__m128 *mass_sse = (__m128*) mass;
	__m128 *xh_sse = (__m128*) xh;
	__m128 *yh_sse = (__m128*) yh;
	__m128 *zh_sse = (__m128*) zh;
	__m128 *xj_sse = (__m128*) xj;
	__m128 *yj_sse = (__m128*) yj;
	__m128 *zj_sse = (__m128*) zj;

	__m128 *axh_sse = (__m128*) axh;
	__m128 *ayh_sse = (__m128*) ayh;
	__m128 *azh_sse = (__m128*) azh;

	float vxj[NPLMAX], vyj[NPLMAX], vzj[NPLMAX], vxh[NPLMAX], vyh[NPLMAX], vzh[NPLMAX];
	float rpl[NPLMAX];

	io_input_planets("pl.in.8", nbod, npl, mass, rpl, xh, yh, zh, vxh, vyh, vzh);
	coord_h2j(nbod, mass, xh, yh, zh, vxh, vyh, vzh, xj, yj, zj, vxj, vyj, vzj);

	test_getacch_sse(nbod, 2, mass_sse, xj_sse, yj_sse, zj_sse, xh_sse, yh_sse, zh_sse);

	io_input_planets("pl.in.8", nbod, npl, mass, rpl, xh, yh, zh, vxh, vyh, vzh);
	coord_h2j(nbod, mass, xh, yh, zh, vxh, vyh, vzh, xj, yj, zj, vxj, vyj, vzj);

	getacch(nbod, mass, xj, yj, zj, xh, yh, zh, axh, ayh, azh);

	//printf("%E\n", mass[1]);
	//printf("%E, %E, %E\n", xh[1], yh[1], zh[1]);
	//printf("%E, %E, %E\n", vxh[1], vyh[1], vzh[1]);
	//printf("%E, %E, %E\n", xj[1], yj[1], zj[1]);
	//printf("%E, %E, %E\n", vxj[1], vyj[1], vzj[1]);
	//printf("%E, %E, %E\n", axh[1], ayh[1], azh[1]);
	//printf("----------\n");
	printf("SCALAR\n");
	printf("%E, %E, %E\n", axh[1], ayh[1], azh[1]);
	printf("%E, %E, %E\n", axh[2], ayh[2], azh[2]);
	printf("%E, %E, %E\n", axh[3], ayh[3], azh[3]);
	printf("%E, %E, %E\n", axh[4], ayh[4], azh[4]);
	printf("%E, %E, %E\n", axh[5], ayh[5], azh[5]);
	printf("%E, %E, %E\n", axh[6], ayh[6], azh[6]);
	printf("%E, %E, %E\n", axh[7], ayh[7], azh[7]);

	printf("-------------------------\n");

	io_input_planets("pl.in.8", nbod, npl, mass, rpl, xh, yh, zh, vxh, vyh, vzh);
	coord_h2j(nbod, mass, xh, yh, zh, vxh, vyh, vzh, xj, yj, zj, vxj, vyj, vzj);

	getacch_sse(2, mass_sse, xj_sse, yj_sse, zj_sse, xh_sse, yh_sse, zh_sse, axh_sse, ayh_sse, azh_sse);

	//printf("%E\n", mass[1]);
	//printf("%E, %E, %E\n", xh[1], yh[1], zh[1]);
	//printf("%E, %E, %E\n", vxh[1], vyh[1], vzh[1]);
	//printf("%E, %E, %E\n", xj[1], yj[1], zj[1]);
	//printf("%E, %E, %E\n", vxj[1], vyj[1], vzj[1]);
	//printf("%E, %E, %E\n", axh[1], ayh[1], azh[1]);
	//printf("----------\n");
	printf("VECTOR\n");
	printf("%E, %E, %E\n", axh[1], ayh[1], azh[1]);
	printf("%E, %E, %E\n", axh[2], ayh[2], azh[2]);
	printf("%E, %E, %E\n", axh[3], ayh[3], azh[3]);
	printf("%E, %E, %E\n", axh[4], ayh[4], azh[4]);
	printf("%E, %E, %E\n", axh[5], ayh[5], azh[5]);
	printf("%E, %E, %E\n", axh[6], ayh[6], azh[6]);
	printf("%E, %E, %E\n", axh[7], ayh[7], azh[7]);

	//io_input_planets("pl.in.8", nbod, npl, mass, rpl, xh, yh, zh, vxh, vyh, vzh);
	//coord_h2j(nbod, mass, xh, yh, zh, vxh, vyh, vzh, xj, yj, zj, vxj, vyj, vzj);

	//double start = omp_get_wtime();
	//for(size_t i = 0; i < 10000000; ++i)
	//{
	//	getacch_sse(2, mass_sse, xj_sse, yj_sse, zj_sse, xh_sse, yh_sse, zh_sse, axh_sse, ayh_sse, azh_sse);
	//}
	//double end = omp_get_wtime();
	//printf("time: %f\n", end - start);

	//io_input_planets("pl.in.8", nbod, npl, mass, rpl, xh, yh, zh, vxh, vyh, vzh);
	//coord_h2j(nbod, mass, xh, yh, zh, vxh, vyh, vzh, xj, yj, zj, vxj, vyj, vzj);

	//start = omp_get_wtime();
	//for(size_t i = 0; i < 10000000; ++i)
	//{
	//	getacch(nbod, mass, xj, yj, zj, xh, yh, zh, axh, ayh, azh);
	//}
	//end = omp_get_wtime();
	//printf("time: %f\n", end - start);


	return 0;


	//coord_h2j(nbod, mass, xh, yh, zh, vxh, vyh, vzh, xj, yj, zj, vxj, vyj, vzj);

	//double start = omp_get_wtime();

	//for(size_t i = 0; i < 10000000; ++i)
	//{
	//	getacch(nbod, mass, xj, yj, zj, xh, yh, zh, axh, ayh, azh);
	//}

	//double end = omp_get_wtime();

	////getacch(nbod, mass, j2rp2, j4rp4, xj, yj, zj, xh, yh, zh, axh, ayh, azh);

	//printf("%E\n", mass[1]);
	//printf("%E, %E, %E\n", xh[1], yh[1], zh[1]);
	//printf("%E, %E, %E\n", vxh[1], vyh[1], vzh[1]);
	//printf("%E, %E, %E\n", xj[1], yj[1], zj[1]);
	//printf("%E, %E, %E\n", vxj[1], vyj[1], vzj[1]);
	//printf("%E, %E, %E\n", axh[1], ayh[1], azh[1]);

	//printf("time: %f\n", end - start);

	//return 0;
}
