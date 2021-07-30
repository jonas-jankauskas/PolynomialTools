**zerocount**
=============

A [FLINT](https://www.flintlib.org/)-based library, written in C,  used for the very fast computation of the number of complex zeros N(f) of polynomials f(x) with rational coefficients (with arbitrarily large numerators and denominators, represented by FLINT fmpq_poly_t type) inside and on the boundary of the unit disk D(0, 1) on the complex plane.

It implements Bistritz algorithm, described in  
Y. Bistritz, *Zero location of polynomials with respect to the unit-circle unhampered by nonessential singularities*, IEEE Transactions on Circuits and Systems I: Fundamental Theory and Applications, **49** (3) (2002), 305-314. See also [Wikipedia article](https://en.wikipedia.org/wiki/Bistritz_stability_criterion).

This library was written and applied to perform the large scale searches for the paper  
K. G. Hare, J. Jankauskas, *On Newman and Littlewood polynomials with a prescribed number of zeros inside the unit disk*, Math. Comp. **90** (2021), 831--870, [arXiv e-print](https://arxiv.org/abs/1910.13994).

The working example can be found in `zerocount.c`

Author:         Jonas Jankauskas  
Copyright:      2021, Jonas Jankauskas. All rights reserved.  
Licence:        LGPLv3  
Version:        1.0  
FLINT version:  2.5.2  
Email:          jonas.jankauskas@gmail.com  
