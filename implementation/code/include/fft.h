#ifndef FFT_H
#define FFT_H

/*
 * Internal functions for Falcon. This is not the API intended to be
 * used by applications; instead, this internal API provides all the
 * primitives on which wrappers build to provide external APIs.
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

#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "fpr.h"

/*
 * "Naming" macro used to apply a consistent prefix over all global
 * symbols.

#ifndef FALCON_PREFIX
#define FALCON_PREFIX   falcon_inner
#endif
#define Zf(name)             Zf_(FALCON_PREFIX, name)
#define Zf_(prefix, name)    Zf__(prefix, name)
#define Zf__(prefix, name)   prefix ## _ ## name  
*/

/*
 * Some computations with floating-point elements, in particular
 * rounding to the nearest integer, rely on operations using _exactly_
 * the precision of IEEE-754 binary64 type (i.e. 52 bits). On 32-bit
 * x86, the 387 FPU may be used (depending on the target OS) and, in
 * that case, may use more precision bits (i.e. 64 bits, for an 80-bit
 * total type length); to prevent miscomputations, we define an explicit
 * function that modifies the precision in the FPU control word.
 *
 * set_fpu_cw() sets the precision to the provided value, and returns
 * the previously set precision; callers are supposed to restore the
 * previous precision on exit. The correct (52-bit) precision is
 * configured with the value "2". On unsupported compilers, or on
 * targets other than 32-bit x86, or when the native 'fpr' type is
 * not used, the set_fpu_cw() function does nothing at all.
 */
static inline unsigned
set_fpu_cw(unsigned x)
{
	return x;
}



/*
 * MSVC 2015 does not know the C99 keyword 'restrict'.
 */
#if defined _MSC_VER && _MSC_VER
#ifndef restrict
#define restrict   __restrict
#endif
#endif

/* ==================================================================== */
/*
 * FFT (falcon-fft.c).
 *
 * A real polynomial is represented as an array of N 'fpr' elements.
 * The FFT representation of a real polynomial contains N/2 complex
 * elements; each is stored as two real numbers, for the real and
 * imaginary parts, respectively. See falcon-fft.c for details on the
 * internal representation.
 */

/*
 * Compute FFT in-place: the source array should contain a real
 * polynomial (N coefficients); its storage area is reused to store
 * the FFT representation of that polynomial (N/2 complex numbers).
 *
 * 'logn' MUST lie between 1 and 10 (inclusive).
 */
void FFT(fpr *f, unsigned logn);

/*
 * Compute the inverse FFT in-place: the source array should contain the
 * FFT representation of a real polynomial (N/2 elements); the resulting
 * real polynomial (N coefficients of type 'fpr') is written over the
 * array.
 *
 * 'logn' MUST lie between 1 and 10 (inclusive).
 */
void iFFT(fpr *f, unsigned logn);

extern const fpr fpr_gm_tab[];

extern const fpr fpr_p2_tab[];

#endif
