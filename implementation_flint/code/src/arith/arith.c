#include "arith.h"
#include "precomputations.h"

/*************************************************
* Name:        arith_setup
*
* Description: Initialize and setup the entire arithmetic backend. 
*              This is strictly required and must be called once
*              before any other function from here is used.
**************************************************/
void arith_setup(void) {
  arith_q_setup();
  arith_real_setup();
  fft_precomp_setup();
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
  arith_real_teardown();
  fft_precomp_teardown();
}

/*************************************************
* Name:        poly_real_from_poly_q
*
* Description: Convert a poly_q into poly_real
*              for arithmetic over R[x]/(x^n + 1) instead
*              of Z[x]/(q, x^n + 1)
* 
* Arguments:   - poly_real res: real-valued polynomial to host the conversion (initialized)
*              - const poly_q arg: polynomial to be converted
**************************************************/
void poly_real_from_poly_q(poly_real res, const poly_q arg) {
  size_t i;
  for (i = 0; i < PARAM_N; i++) {
    poly_real_set_si(
      res,
      i,
      poly_q_get_coeff_centered(arg, i)
    );
  }
}

/*************************************************
* Name:        poly_q_from_poly_real
*
* Description: Convert a poly_real into poly_q by rounding
*              for arithmetic over Z[x]/(q, x^n + 1) instead
*              of R[x]/(x^n + 1)
* 
* Arguments:   - poly_q res: polynomial to host the conversion (initialized)
*              - const poly_real arg: real-valued polynomial to be rounded/converted
**************************************************/
void poly_q_from_poly_real(poly_q res, const poly_real arg) {
  size_t i;
  for (i = 0; i < PARAM_N; i++) {
    poly_q_set_coeff(
      res,
      i,
      poly_real_get_coeff_rounded(arg, i)
    );
  }
}

/*************************************************
* Name:        poly_q_samplefz
*
* Description: Gaussian sampler SampleFz over Z[x]/(x^n + 1) from
*              [GM18] and conversion to poly_q
* 
* Arguments:   - poly_q res: polynomial to host the Gaussian sample (initialized)
*              - const poly_real f: variance polynomial for the Gaussian (Careful: not the stddev)
*              - const poly_real c: center polynomial for the Gaussian
**************************************************/
void poly_q_samplefz(poly_q res, const poly_real f, const poly_real c) {
  poly_real tmp;
  poly_real_init(tmp);
  poly_real_samplefz(tmp, f, c);
  poly_q_from_poly_real(res, tmp);
  poly_real_clear(tmp);
}
