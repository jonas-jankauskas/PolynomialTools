/** @file		zerocount.c
	@brief		source code for Bistritz_rule example
	@author		Jonas Jankauskas
	@date		November 22, 2019
	@version	1.0
	@note		Flint version 2.5.2
*/

#include "flint/flint.h"
#include "flint/fmpq_poly.h"

#include "debug.h"
#include "zerocount.h"

/**
	@brief		Main function
	@param		argc - int, number of command line parameters 
	@param		argv[] - char *, pointer to null-terminated strings holding the command line params
	@return		nothing
*/

int main(int argc, char *argv[]) {
	
	//debug routine
	DEBUG_ENTER_AT(2);
	
	//FLINT rational polynomial
	fmpq_poly_t P;
	
	//root numbers
	slong roots_iuc, roots_uc;
	
	// initialize polynomial and read its coefficients
	fmpq_poly_init(P);
	flint_printf ("# P(x):");
	fmpq_poly_read(P);
	
	// calculate root numbers
	Bistritz_rule(&roots_iuc, &roots_uc, P);
	
	// print the answer
	flint_printf("#Got P(x)= ");
	fmpq_poly_print_pretty(P, "x");
	flint_printf("\n# Zeros inside/on the unit circle:\n%wd %wd\n", roots_iuc, roots_uc);
	
	//clear allocated memory
	fmpq_poly_clear(P);
	
	//debug routine
	DEBUG_LEAVE_AT(2);
	
}