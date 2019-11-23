/** @file		bistritz.c
	@brief		source code for Bistritz_rule and auxiliary functions
	@author		Jonas Jankauskas
	@date		November 22, 2019
	@version	1.0
	@note		FLINT version 2.5.2
	@details	Contains the source code of the Bistritz algorithm, together with routines that perform calculations symmetric polynomials (saving about half of computational time) and 3-term Bistritz recurrence formula. Reasonable attempt to optimize the code for the speed was made without compromising the code readability. The program is implemented via FLINT. FLINT conventions are followed where possible: each function returns its result via the 1st argument. It is assumed that function arguments are initialized by user by appropriate FLINT initialization routines before the call. Unless the variable is declared const, users must assume that the function will modify the contents of that variable.
*/

#include "flint/flint.h"
#include "flint/fmpz.h"
#include "flint/fmpq.h"
#include "flint/fmpq_poly.h"

#include "debug.h"


/**
	@brief		Returns the index of lowest non-zero term;
	@details	0 if poly = 0
	@param		poly FLINT fmpq_poly_t type rational polynomial
	@return		FLINT slong type that is equal to index of first non-zero term in poly	
*/
slong get_lambda(const fmpq_poly_t poly) {
	
	/* debug */
	DEBUG_ENTER_AT(3);
	
	for (slong i=0; i < poly->length; i++)
		if (!fmpz_is_zero(poly->coeffs+i)) {
			
			/* debug */
			DEBUG_LEAVE_AT(3);
			
			return(i);
		}
	
	/* debug */
	DEBUG_LEAVE_AT(3);
	
	return(0);
	
}

/**
	@brief	Returns reference to numerator of n-th coefficient or 0 if index n is out of range
	@param	poly - FLINT fmpq_poly_t type rational polynomial
	@param	n - FLINT slong type index
	@param	const_null - a reference to the FLINT fmpz_t variable that contains 0 constant
	@return	reference to the numerator of n-th coefficient of poly as a pointer to FLINT fmpz	
*/
const fmpz *get_coeff_numref(const fmpq_poly_t poly,
						const fmpz_t const_null,
						const slong n) {
							
	/* debug */
	DEBUG_ENTER_AT(4);
	
	if ((n >= 0) && (n < poly->length)) {
		
		/* debug */
		DEBUG_LEAVE_AT(4);
	
		return(poly->coeffs+n);
		
	}
	else {
		/* debug */
		DEBUG_LEAVE_AT(4);
		
		return(const_null);
	
	}
	
}

/**
	@brief		Returns n-th coefficient of poly or 0 if index n is out of range
	@details	Same as fmpq_poly_get_coeff_fmpq with a check for negative n
	@param		coeff - FLINT fmpq_t type rational coefficient
	@param		poly - FLINT fmpq_poly_t type rational polynomial
	@param		n - FLINT slong type index
	@return		result is returned in coeff	
*/
void get_coeff(fmpq_t coeff, const fmpq_poly_t poly, const slong n) {
	
	/* debug */
	DEBUG_ENTER_AT(4);
	
	if ((n >= 0) && (n < poly->length)) {
		fmpq_set_fmpz_frac(coeff, poly->coeffs+n, poly->den);
		fmpq_canonicalise(coeff);
	}
	else
		fmpq_zero(coeff);
	
	/* debug */
	DEBUG_LEAVE_AT(4);
	
	return;
}

/**
	@brief		Calculates delta coefficient for the three-term recursion
	@param		delta - recursion coefficient, FLINT fmpq_t type 
	@param		T_prev - previous polynomial in the recursion, FLINT fmpq_poly_t type.
	@param		T_curr - current polynomial in the recursion, FLINT fmpq_poly_t type.
	@return		result is returned in delta	
*/
void get_delta(fmpq_t delta,
		const fmpq_poly_t T_prev,
		const fmpq_poly_t T_curr) {
			
	/* debug */
	DEBUG_ENTER_AT(3);
	
	if ((fmpq_poly_length(T_prev) > 0) && (fmpq_poly_length(T_curr) > 0)) {
		
		get_coeff(delta, T_curr, get_lambda(T_curr));
		fmpz_swap(&(delta->num), &(delta->den));
		fmpz_mul(&(delta->num), &(delta->num), T_prev->coeffs);
		fmpz_mul(&(delta->den), &(delta->den), T_prev->den);
		fmpq_canonicalise(delta);
	}
	else
		fmpq_zero(delta);
	
	DEBUG_LEAVE_AT(3);
	
	return;
} 

/**
	@brief		Divides antisymmetric P(x) (P*=-P) by (x-1) in  Q[x], returns quotient
	@details	Horner evaluation at 1 and exploitation of negative-symmetry
	@param		res - quotient polynomial, FLINT fmpq_poly_t type
	@param		poly - rational polynomial that is to be divided, FLINT fmpq_poly_t type
	@return		result is returned in res	
*/
void div_x_minus_1_asym(fmpq_poly_t res, const fmpq_poly_t poly) {
	
	slong deg, length, start, mid;
	
	/* debug */
	DEBUG_ENTER_AT(3);
	
	fmpq_poly_set(res, poly);
	
	deg = fmpq_poly_degree(res);
	
	start =  get_lambda(res);
	length = fmpq_poly_length(res)- start;
	mid = length/2+start;
	
	/* Since P=-P*, the quotient P/(x-1) is symmetric */
		
	for (slong i = deg - 1; i >= mid; i--)	//calculate first half of the result and the middle term
			fmpz_add(res->coeffs + i, res->coeffs + i, res->coeffs + (i+1));
		
	for (slong i = 1; i < length/2; i++) //copy the second half using the symmetry
			fmpz_set(res->coeffs + start+i, res->coeffs + (deg+1-i));
			
	if (start >= 1)
		fmpz_zero(res->coeffs+start);
		
	fmpq_poly_shift_right(res, res, 1);	//discard last term
	fmpq_poly_canonicalise(res);
	
	/* debug */
	DEBUG_LEAVE_AT(3);
	
	return;

}

/* ------ Evaluates formaly-symmetric polynomial at x=1 ------ */
/**
	@brief		Evaluates formaly-symmetric polynomial at x=1
	@details	Uses Horner evaluation at 1 and symmetry of P
	@param		val - the computed FLINT fmpq_t value P(1)
	@param		poly - rational polynomial of FLINT fmpq_poly_t type
	@return		result is returned in val	
*/

void eval_at_1_sym(fmpq_t val, const fmpq_poly_t poly) {
	
	slong length, start, mid;
	
	/* debug */
	DEBUG_ENTER_AT(3);
	
	fmpq_zero(val);
	
	if (fmpq_poly_length(poly) > 0) {
		
		start = 0;
		while (fmpz_is_zero(poly->coeffs+start)) start++; //first non-zero term
		
		length = fmpq_poly_length(poly) - start;
		mid = start+length/2;
	
		for (slong i = start; i < mid; i++)				//sums the first half
			fmpz_add(&(val->num), &(val->num), poly->coeffs+i);
		
		fmpz_mul_2exp(&(val->num), &(val->num), 1);		//doubles the sum
		
		if ((length % 2) != 0)							//account for the middle term
			fmpz_add(&(val->num), &(val->num), poly->coeffs+mid);
		
		fmpz_set(&(val->den), poly->den);
		fmpq_canonicalise(val);
	}
	
	return;
	
	/* debug */
	DEBUG_LEAVE_AT(3);
	
	}

/**
	@brief		Repeatedly divides P(x) by (x-1) until P(-1)<>0, returns quotient.
	@details	Quotient polynomial is obtained by a Horner-like evaluation at x=1
	@param		count - multiplicity of root x=1 in P(x). If P(x)=0, 0 is returned.
	@param		poly - FLINT fmpq_poly_t type rational polynomial
	@return		result is returned in count	
*/

	
void clear_x_minus_1(slong *count, fmpq_poly_t poly) {
	
	fmpq_poly_t P;
	
	int cont;
	
	/* debug */
	DEBUG_ENTER_AT(2);
	
	*count = 0;
	
	if (fmpq_poly_length(poly) > 0) {
		
		fmpq_poly_init(P);
		fmpq_poly_set(P, poly);
		
		cont = 1;
		
		do {
			//Quotient polynomal by Horner-like evaluation at 1;
			for (slong i = fmpq_poly_degree(P)-1; i >= 0; i--)
				fmpz_add(P->coeffs+i, P->coeffs+i, P->coeffs+i+1);
			
			if (fmpz_is_zero(P->coeffs)) {
				
				(*count)++;
				fmpq_poly_shift_right(P, P, 1); //divide & forget the remainder
				fmpq_poly_canonicalise(P);
				fmpq_poly_set(poly, P);			//the result so far
			}
			else cont = 0;
			
		} while (cont);
		fmpq_poly_clear(P);
	}
	
	return;
	
	/* debug */
	DEBUG_LEAVE_AT(2);
	}

/**
	@brief		Initializes polynomials T1, T2 and their signs at x=1
	@details	T1 = (D+D*), T2 = (D-D*)/(z-1), sigma1 = sgn(T1(1)), sigma2 = sgn(T2(1))
	@param		T1 - first polynomial of the recurrence, FLINT fmpq_poly_t type
	@param		T2 - second polynomial of the recurrence, FLINT fmpq_poly_t type
	@param		sigma1 - sign of T1 at x=1, fmpq_t type rational number
	@param		sigma2 - sign of T2 at x=1, fmpq_t type rational number
	@param		D - initial FLINT fmpq_poly_t type rational polynomial
	@return		results are is stored T1, T2, sigma1, sigma2	
*/

void rule_init(fmpq_poly_t T1, fmpq_poly_t T2, fmpq_t sigma1, fmpq_t sigma2, const fmpq_poly_t D) {
	
	/* (re)-initialization */
	
	/* debug */
	DEBUG_ENTER_AT(2);
	
	fmpq_poly_reverse(T1, D, fmpq_poly_length(D));
	
	fmpq_poly_sub(T2, D, T1);					
	div_x_minus_1_asym(T2, T2);
	
	fmpq_poly_add(T1, T1, D);
	
	eval_at_1_sym(sigma1, T1);
	eval_at_1_sym(sigma2, T2);
	
	/* debug */
	DEBUG_LEAVE_AT(2);
	
	return;
	
}

/**
	@brief		Debuging function; prints the recursion data to the stdout
	@details	Prints n-th Polynomial T_n, lambda value of T_n and sigma=sgn(Tn(1))
	@param		n - FLINT slong index of a sequence member
	@param		T - n-th polynomial of the recurrence, FLINT fmpq_poly_t type
	@param		sigma - sign of T_n at x=1, fmpq_t type rational number	
*/
void print_T_data(slong n, fmpq_poly_t T, fmpq_t sigma) {
	
	/* debug */
	DEBUG_ENTER_AT(4);
	
	flint_printf("# T_%wd = ", n);
	fmpq_poly_print(T);
	flint_printf("\n");
	
	flint_printf("# lambda_%wd = %wd, sigma_%wd = ", n, get_lambda(T), n);
	fmpq_print(sigma);
	flint_printf("\n");
	
	/* debug */
	DEBUG_LEAVE_AT(4);
	
	return;
	
}

/* ------ Singular case: if T2 = 0, T1 <> 0, then re-initialize T2 and T3 from T1' ------ */

/**
	@brief		Re-initializes polynomials T2, T3 and their signs at x=1 in singular cases
	@details	D = T1', T2 = (D+D*), T3 = -(D-D*)/(z-1), sigma1 = sgn(T1(1)), sigma2 = sgn(T2(1)), sigma3 = sgn(T3(1))
	@param		T1 - first polynomial of the recurrence, FLINT fmpq_poly_t type; used as aplaceholder for T1'
	@param		T2 - second polynomial of the recurrence, FLINT fmpq_poly_t type
	@param		T3 - third polynomial of the recurrence, FLINT fmpq_poly_t type
	@param		sigma1 - sign of T1 at x=1, fmpq_t type rational number
	@param		sigma2 - sign of T2 at x=1, fmpq_t type rational number
	@param		sigma3 - sign of T3 at x=1, fmpq_t type rational number
	@return		results are stored in T1, T2, T3, sigma1, sigma2, sigma3	
*/

void do_singular(fmpq_poly_t T1,
				 fmpq_poly_t T2,
				 fmpq_poly_t T3,
				 fmpq_t sigma1,
				 fmpq_t sigma2,
				 fmpq_t sigma3) {
	
	/* debug */
	DEBUG_ENTER_AT(2);
	
	// Singular case: re-initialize from T'(z) of  a last nonzero T(z) <> 0
	fmpq_poly_derivative(T1, T1);
	// As '*' differs from Bistritz '#' operation: initialize from the T'(z), then negate
	rule_init(T2, T3, sigma2, sigma3, T1);
	fmpq_poly_neg(T2, T2);
	fmpq_neg(sigma2, sigma2);
	
	/* debug */
	DEBUG_LEAVE_AT(2);
	
	return;
	
}

/**
	@brief		Calculates next polynomial T3 in regular recursion
	@details	T3 = delta*(z^(2*lambda+1)+1)*T_2*z^(-lambda) - T1
	@param		T3 - new polynomial of the recurrence, FLINT fmpq_poly_t type
	@param		sigma3 - sign of next polynomial T3 at x=1, fmpq_t type rational number
	@param		delta - placeholder for the recursion coefficient
	@param 		fmpz_const_null - fmpz_t placeholder for '0'
	@param		T1 - pre-previous polynomial of the recurrence, FLINT fmpq_poly_t type
	@param		T2 - previous polynomial of the recurrence, FLINT fmpq_poly_t type
	@param		sigma1 - sign of T1 at x=1, fmpq_t type rational number
	@param		sigma2 - sign of T2 at x=1, fmpq_t type rational number
	@param		flength - FLINT slong (formal) symmetrical length of the next polynomial, used in #-operation
	@return		results are stored in T3, sigma3, delta	
*/

void do_recurence(fmpq_poly_t T3, 					//next polynomial T
					fmpq_t sigma3,					//next sigma
					fmpq_t delta, 					//delta coefficient and temporary variable
					const fmpz_t fmpz_const_null,	//placeholder for '0'
					const fmpq_poly_t T1,
					const fmpq_poly_t T2,
					const fmpq_t sigma1,
					const fmpq_t sigma2,
					const slong flength) {			//(formal) symmetrical length of the next polynomial
		
	slong mid, lambda;
	
	/* debug */
	DEBUG_ENTER_AT(2);
	
	//Handle 0-valued polynomials separately, because FLINT threats them as empty 
	if (flength <= 0) {
		
		fmpq_poly_zero(T3);
		fmpq_zero(sigma3);
		
		/* debug */
		DEBUG_LEAVE_AT(2);
		
		return;
	}
	
	//coefficients
	lambda = get_lambda(T2);
	get_delta(delta, T1, T2);
	
	/* Debug */
	DEBUG_MSG_AT(2, "# regular case: flength=%wd, lambda=%wd, ", flength, lambda);
	DEBUG_MSG_AT(2, "delta=");
	DEBUG_FMPQ_AT(2, delta);
	DEBUG_MSG_AT(2, "\n");
	
	//next sigma
	fmpq_mul(sigma3, delta, sigma2);
	fmpq_mul_2exp(sigma3, sigma3, 1);
	fmpq_sub(sigma3, sigma3, sigma1);
	
	//enough place for T3?
	if (T3->alloc < flength) {			
		fmpq_poly_realloc(T3, flength);
	}
	T3->length = flength;
	
	//T3 middle term
	mid = (flength-1)/2;
	
	/*	Debug*/
	DEBUG_MSG_AT(2, "# allocated=%wd length=%wd mid=%wd\n", T3->alloc, T3->length, mid);
	
	// three-term recursion 
	// t3[i] = delta*(t2[i-lambda] + t2[i+lambda+1]) - t1[i+1]
	
	//pre-multiply numerators and denominators - delta is used as a placeholder!
	fmpz_mul(&(delta->num), &(delta->num), T1->den);
	fmpz_mul(&(delta->den), &(delta->den), T2->den);
	
	for (slong i=0; i <= mid; i++) {
		
		fmpz_add(T3->coeffs+i, get_coeff_numref(T2, fmpz_const_null, i-lambda), get_coeff_numref(T2, fmpz_const_null, i+lambda+1));
		fmpz_mul(T3->coeffs+i, T3->coeffs+i, &(delta->num));
		fmpz_submul(T3->coeffs+i, &(delta->den), get_coeff_numref(T1, fmpz_const_null, i+1));
		
	}
	
	//by symmetry of T3
	for (slong i=mid+1; i < flength; i++)
		fmpz_set(T3->coeffs+i, T3->coeffs+flength-i-1);
	//possible latter optimization?
	
	//denominator of T3
	fmpz_mul(T3->den, &(delta->den), T1->den);
	fmpq_poly_canonicalise(T3);
	
	/* debug */
	DEBUG_LEAVE_AT(2);
	
	return;
	
}

/**
	@brief		Zero counting procedure
	@details	If poly <> 0, returns no. of complex zeros of poly inside |z| < 1 and on |z| = 1 through *in_uc and *on_uc. If poly == 0, stores the pair of numbers -1, 0.
	@param		poly - FLINT fmpq_poly_t type rational polynomial (must be initialized before the call)
	@param		in_uc - pointer to FLINT slong, no. of complex zeros inside unit circle |z|<1
	@param		on_uc - pointer to FLINT slong, no. of complex zeros on unit circle |z|=1
	@return		results are returned through in_uc, on_uc variables
*/ 
void Bistritz_rule(slong *in_uc, slong *on_uc, const fmpq_poly_t poly) {
	
	fmpz_t const_null;
	
	fmpq_t sigma_prev, sigma_curr, sigma_next, delta;
	
	fmpq_poly_t D, T_prev, T_curr, T_next;
	
	slong deg, last_sgn, curr_sgn, singular, vars, vars_reg;
	
	/* debug */
	DEBUG_ENTER_AT(1);
	
	fmpz_init(const_null);
	fmpz_zero(const_null);
	
	fmpq_poly_init(D);
	fmpq_poly_init(T_prev);
	fmpq_poly_init(T_curr);
	fmpq_poly_init(T_next);
	
	fmpq_init(delta);
	fmpq_init(sigma_prev);
	fmpq_init(sigma_curr);
	fmpq_init(sigma_next);
	
	vars = 0;
	vars_reg = 0;
	singular = -1;
	
	/* debug */
	DEBUG_MSG_AT(1, "# received:\n# poly = ");
	DEBUG_FMPQ_POLY_AT(1, poly);
	DEBUG_MSG_AT(1, "\n");
	
	fmpq_poly_set(D, poly);
	clear_x_minus_1(on_uc, D);
	
	deg = fmpq_poly_degree(D);
	
	/* debug */
	DEBUG_MSG_AT(1, "# (x-1) factors cleared, degree deg=%wd\n# D = ", deg);
	DEBUG_FMPQ_POLY_AT(1, poly);
	DEBUG_MSG_AT(1, "\n");
	
	rule_init(T_prev, T_curr, sigma_prev, sigma_curr, D);
	last_sgn = fmpq_sgn(sigma_prev);
	
	/* debug */
	DEBUG_T(1, deg, T_prev, sigma_prev);
	
	for (slong i = deg-1; i >= 0; i--) {
		
		/* debug */
		DEBUG_MSG_AT(1, "# * loop i = %wd *\n", i);
		DEBUG_T(1, i, T_curr, sigma_curr);
		
		if (fmpq_poly_is_zero(T_curr)) {
			
			if (fmpq_poly_is_zero(T_prev))
				break;	
			
			else {
				if (fmpz_is_zero(T_prev->coeffs))
					do_recurence(T_next, sigma_next, delta, const_null, T_prev, T_curr, sigma_prev, sigma_curr, i);
				else { 
			
					// in singular case: re-initialize the recursion from T_{k+1}'(z) <> 0
					do_singular(T_prev, T_curr, T_next, sigma_prev, sigma_curr, sigma_next);
					
					/* debug */
					DEBUG_MSG_AT(1, "# singularity after s=%wd:\n", i+1);
					DEBUG_T(1, i, T_curr, sigma_curr);
			
					//record position of the first singularity and sign variations before it
					if (singular == -1) {
						singular = i;
						vars_reg = vars;
						
						DEBUG_MSG_AT(1, "# vars_reg = %wd sign variations occured before singularity.\n", vars_reg);
					}
				}
			}
		}
		else
			do_recurence(T_next, sigma_next, delta, const_null, T_prev, T_curr, sigma_prev, sigma_curr, i);
		
		/* Count sign variations */
		curr_sgn = fmpq_sgn(sigma_curr);
		vars += (1-last_sgn*curr_sgn)/2;
		if (curr_sgn != 0)
			last_sgn = curr_sgn;
		
		/* Step down from i to i-1 */
		fmpq_poly_swap(T_prev, T_curr);
		fmpq_poly_swap(T_curr, T_next);
		fmpq_swap(sigma_prev, sigma_curr);
		fmpq_swap(sigma_curr, sigma_next);

	}
	
	/* debug */
	DEBUG_MSG_AT(1, "# * end loop *\n");
	DEBUG_T(1, -1, T_curr, sigma_curr); //last produced polynomial
	
	if (singular== -1)
		vars_reg = vars;
	
	*in_uc = deg - vars;
	*on_uc += 2*(vars-vars_reg) - singular - 1;
	
	/* debug */
	DEBUG_MSG_AT(1, "# singular =%wd/vars_reg=%wd/vars=%wd\n",  singular, vars_reg, vars);
	DEBUG_MSG_AT(1, "# roots IUC/UC: %wd/%wd\n", *in_uc, *on_uc); //last produced polynomial
	
	fmpq_poly_clear(D);
	fmpq_poly_clear(T_prev);
	fmpq_poly_clear(T_curr);
	fmpq_poly_clear(T_next);
	
	fmpq_clear(sigma_prev);
	fmpq_clear(sigma_curr);
	fmpq_clear(sigma_next);
	fmpq_clear (delta);
	
	fmpz_clear(const_null);
	
	DEBUG_LEAVE_AT(1);
	
	return;
}