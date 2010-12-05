/*
 * corrector.h
 *
 *  Created on: Nov 23, 2010
 *      Author: Huey
 */

#ifndef CORRECTOR_H_
#define CORRECTOR_H_

void invcorrector(NBodies bodies, float tinc, float dt);
void correctterm_hel(NBodies bodies);
void correctterm_bar(NBodies bodies);


#endif /* CORRECTOR_H_ */
