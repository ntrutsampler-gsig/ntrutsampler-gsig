#include "arith.h"
#include "fft.h"
#include "fpr.h"

/*************************************************
* Name:        arith_setup
*
* Description: Initialize and setup the entire arithmetic backend. 
*              This is strictly required and must be called once
*              before any other function from here is used.
**************************************************/
void arith_setup(void) {
  arith_q_setup();
}

/*************************************************
* Name:        arith_teardown
*
* Description: Clean up and teardown the entire arithmetic backend. 
*              This is strictly required and must be called once at 
*              the very end to release any resources.
**************************************************/
void arith_teardown(void) {
  arith_q_teardown();
}

/*************************************************
* Name:        poly_real_from_poly_q
*
* Description: Convert a poly_q into poly_real in FFT representation
*              for arithmetic over R[x]/(x^n + 1) instead
*              of Z[x]/(q, x^n + 1) 
* 
* Arguments:   - poly_real res: real-valued polynomial to host the conversion (initialized)
*              - const poly_q arg: polynomial to be converted
**************************************************/
void poly_real_from_poly_q(poly_real res, const poly_q arg) {
  size_t i;
  for (i = 0; i < PARAM_N; i++) {
    poly_real_set_si(res, i, poly_q_get_coeff_centered(arg, i));
  }
  FFT(res, PARAM_LOG_N);
  // WARNING: Results in polynomial in FFT representation
}

/*************************************************
* Name:        poly_q_from_poly_real
*
* Description: Convert a poly_real in FFT representation into poly_q by rounding
*              for arithmetic over Z[x]/(q, x^n + 1) instead
*              of R[x]/(x^n + 1)
* 
* Arguments:   - poly_q res: polynomial to host the conversion (initialized)
*              - const poly_real arg: real-valued polynomial to be rounded/converted
**************************************************/
void poly_q_from_poly_real(poly_q res, const poly_real arg) {
  size_t i;
  // WARNING: Requires polynomial in FFT representation
  fpr tmp[PARAM_N];
  for (i = 0; i < PARAM_N; i++) {
    tmp[i] = arg[i];
  }
  iFFT(tmp, PARAM_LOG_N);
  for (i = 0; i < PARAM_N; i++) {
    poly_q_set_coeff(res, i, poly_real_get_coeff_rounded(tmp, i));
  }
}