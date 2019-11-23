/** @file		zerocount.h
	@author		Jonas Jankauskas
	@date		November 22, 2019
	@version	1.0
	@note		Flint version 2.5.2
	@brief		Header file for the polynomial zero counting function prototypes
	@details	Currently contains only the prototype of the Bistritz algorithm procedure, used to count the number of complex zeros of a polynomial P(x) of Q[x] inside and on the boundary of complex unit circle Disk(0, 1) = {z: |z| <= 1}.
*/

#ifndef ZEROCOUNT_H
	
	#define ZEROCOUNT_H
	#include "flint/fmpq_poly.h"

	void Bistritz_rule(slong *, slong *, const fmpq_poly_t);


#endif