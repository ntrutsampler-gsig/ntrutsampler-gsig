#ifndef FPR_H
#define FPR_H

#include <stdint.h>
#include <stdlib.h>
#include "params.h"

/*
 * Floating-point operations.
 *
 * ==========================(LICENSE BEGIN)============================
 *
 * Copyright (c) 2017-2019  Falcon Project
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 * ===========================(LICENSE END)=============================
 *
 * @author   Thomas Pornin <thomas.pornin@nccgroup.com>
 */

#include <math.h>

typedef double fpr;

/* ==================================================================== */
/*
 * Wrapper of floating-point real numbers (fpr.h, fpr.c) using C double.
 */

/*
 * Real numbers are implemented by an extra header file, included below.
 * This is meant to support pluggable implementations. The default
 * implementation relies on the C type 'double'.
 *
 * The included file must define the following types, functions and
 * constants:
 *
 *   fpr
 *         type for a real number
 *
 *   fpr fpr_of(int64_t i)
 *         cast an integer into a real number; source must be in the
 *         -(2^63-1)..+(2^63-1) range
 *
 *	 int64_t fpr_floor(fpr x)
 *		   compute floor(x)
 *
 *   fpr fpr_add(fpr x, fpr y)
 *         compute x + y
 *
 *   fpr fpr_sub(fpr x, fpr y)
 *         compute x - y
 *
 *   fpr fpr_neg(fpr x)
 *         compute -x
 *
 *   fpr fpr_half(fpr x)
 *         compute x/2
 *
 *   fpr fpr_double(fpr x)
 *         compute x*2
 *
 *   fpr fpr_mul(fpr x, fpr y)
 *         compute x * y
 *
 *   fpr fpr_sqr(fpr x)
 *         compute x * x
 *
 *   fpr fpr_div(fpr x, fpr y)
 *         compute x/y
 *
 *   fpr fpr_inv(fpr x)
 *         compute 1/x
 *
 *   fpr fpr_sqrt(fpr x)
 *         compute the square root of x
 *
 *   int fpr_lt(fpr x, fpr y)
 *         return 1 if x < y, 0 otherwise
 *
 *	 fpr fpr_expm(fpr x)
 *		   return exp(-x)
 *
 *	 int64_t fpr_ceil(fpr x)
 *		   compute ceil(x)
 *
 *	 int fpr_leq(fpr x, fpr y)
 *		   return 1 if x <= y, 0 otherwise
 *
 *   const fpr fpr_gm_tab[]
 *         array of constants for FFT / iFFT
 *
 *   const fpr fpr_p2_tab[]
 *         precomputed inverses of half powers of 2 (2/2^i) (by index, i = 0 to 10)
 */

static inline fpr
fpr_of(int64_t i)
{
	return (fpr)i;
}

static inline int64_t
fpr_floor(fpr x)
{
	return (int64_t)floor(x);
}

static inline fpr 
fpr_add(fpr x, fpr y)
{
	return x + y;
}

static inline fpr
fpr_sub(fpr x, fpr y)
{
	return x - y;
}

static inline fpr
fpr_neg(fpr x)
{
	return -x;
}

static inline fpr
fpr_half(fpr x)
{
	return x / (double)2;
}

static inline fpr
fpr_double(fpr x)
{
	return 2 * x;
}

static inline fpr 
fpr_mul(fpr x, fpr y)
{
	return x * y;
}

static inline fpr
fpr_sqr(fpr x)
{
	return fpr_mul(x, x);
}

static inline fpr 
fpr_div(fpr x, fpr y)
{
	return x / y;
}

static inline fpr
fpr_inv(fpr x)
{
	return fpr_div(PARAM_ONE, x);
}

static inline fpr 
fpr_sqrt(fpr x)
{
	return sqrt(x);
}

static inline int
fpr_lt(fpr x, fpr y)
{
	return (int)(x < y);
}

static inline fpr 
fpr_expm(fpr x)
{
	return exp(-x);
}

static inline int64_t
fpr_ceil(fpr x)
{
	return (int64_t)ceil(x);
}

static inline int64_t
fpr_leq(fpr x, fpr y)
{
	/*
	 * We use that [x <= y] = NOT([y < x]) = ([y < x]&1) ^ 1
	 */
	return (fpr_lt(y,x) & 1) ^ 1;
}

/* ====================================================================== */

extern const fpr fpr_gm_tab[];

extern const fpr fpr_p2_tab[];

#endif
