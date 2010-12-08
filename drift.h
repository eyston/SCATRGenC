/*
 * drift.h
 *
 *  Created on: Nov 23, 2010
 *      Author: Huey
 */

#ifndef DRIFT_H_
#define DRIFT_H_

#include "structures.h"

void drift(NBodies bodies, float dt);
void drift_one(float mu, float xj, float yj, float zj, float vxj, float vyj, float vzj, float dt);
void drift_dan(float mu, float xj, float yj, float zj, float vxj, float vyj, float vzj, float dt);
void drift_kepmd(float dm, float es, float ec, float xkep, float s, float c);
void drift_kepu(float dt, float r0, float mu, float alpha, float u, float rp, float c1, float c2, float c3, int iflg);
void drift_kepu_guess(float dt, float r0, float mu, float alpha, float u, float s);
void drift_kepu_new(float s, float dt, float r0, float mu, float alpha, float u, float fp, float c1, float c2, float c3, int iflg);
void drift_kepu_fchk(float dt, float r0, float mu, float alpha, float u, float st, float fo);
void drift_kepu_lag(float s, float dt, float r0, float mu, float alpha, float u, float fp, float c1, float c2, float c3, float iflg);
void drift_kepu_p3solve(float dt, float r0, float mu, float alpha, float u, float s, int iflg);
void drift_kepu_stumpff(float x, float c0, float c1, float c2, float c3);

#endif /* DRIFT_H_ */
