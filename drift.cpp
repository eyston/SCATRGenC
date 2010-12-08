/*
 * drift.cpp
 *
 *  Created on: Nov 23, 2010
 *      Author: Huey
 */

#include <stdio.h>
#include "structures.h"
#include "drift.h"

void drift(NBodies bodies, float dt)
{
	printf("drift called\n");
	for(int i = 1; i < bodies.nbod; i++)
	{
		float mu = 0.0f;
		drift_one(mu, bodies.xj[i], bodies.yj[i], bodies.zj[i], bodies.vxj[i], bodies.vyj[i], bodies.vzj[i], dt);
	}
}

void drift_one(float mu, float xj, float yj, float zj, float vxj, float vyj, float vzj, float dt)
{
	printf("drift_one called\n");

	drift_dan(mu, xj, yj, zj, vxj, vyj, vzj, dt);
}

void drift_dan(float mu, float xj, float yj, float zj, float vxj, float vyj, float vzj, float dt)
{
	printf("drift_dan called\n");
	float dm = 0.0f, es = 0.0f, ec = 0.0f, xkep = 0.0f, s = 0.0f, c = 0.0f;
	drift_kepmd(dm, es, ec, xkep, s, c);
	float r0 = 0.0f, alpha = 0.0f, u = 0.0f, rp = 0.0f, c1 = 0.0f, c2 = 0.0f, c3 = 0.0f;
	int iflg = 0;
	drift_kepu(dt, r0, mu, alpha, u, rp, c1, c2, c3, iflg);
}

void drift_kepmd(float dm, float es, float ec, float xkep, float s, float c)
{
	// leaf function
	printf("drift_kepmd called\n");
}

void drift_kepu(float dt, float r0, float mu, float alpha, float u, float rp, float c1, float c2, float c3, int iflg)
{
	printf("drift_kepu called\n");

	float s = 0.0f;

	drift_kepu_guess(dt, r0, mu, alpha, u, s);

	float fp = 0.0f;

	drift_kepu_new(s, dt, r0, mu, alpha, u, fp, c1, c2, c3, iflg);

	float fo = 0.0f, st = 0.0f;

	drift_kepu_fchk(dt, r0, mu, alpha, u, st, fo);

	drift_kepu_lag(s, dt, r0, mu, alpha, u, fp, c1, c2, c3, iflg);

}

void drift_kepu_guess(float dt, float r0, float mu, float alpha, float u, float s)
{
	printf("drift_kepu_guess called\n");

	int iflg = 0;

	drift_kepu_p3solve(dt, r0, mu, alpha, u, s, iflg);
}

void drift_kepu_new(float s, float dt, float r0, float mu, float alpha, float u, float fp, float c1, float c2, float c3, int iflg)
{
	printf("drift_kepu_new called\n");

	float x = 0.0f, c0 = 0.0f;

	drift_kepu_stumpff(x, c0, c1, c2, c3);
}

void drift_kepu_fchk(float dt, float r0, float mu, float alpha, float u, float st, float fo)
{
	printf("drift_kepu_fchk called\n");

	float x = 0.0f, c0 = 0.0f, c1 = 0.0f, c2 = 0.0f, c3 = 0.0f;

	drift_kepu_stumpff(x, c0, c1, c2, c3);
}

void drift_kepu_lag(float s, float dt, float r0, float mu, float alpha, float u, float fp, float c1, float c2, float c3, float iflg)
{
	printf("drift_kepu_lag called\n");

	float x = 0.0f, c0 = 0.0f;

	drift_kepu_stumpff(x, c0, c1, c2, c3);
}

void drift_kepu_p3solve(float dt, float r0, float mu, float alpha, float u, float s, int iflg)
{
	// leaf function
	printf("drift_kepu_p3solve called\n");
}

void drift_kepu_stumpff(float x, float c0, float c1, float c2, float c3)
{
	// leaf function
	printf("drift_kepu_stumpff called\n");
}
