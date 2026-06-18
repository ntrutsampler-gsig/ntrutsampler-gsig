#include "poly_real.h"
#include "random.h"
#include "macros.h"
#include "fft.h"
#include <math.h>

/*************************************************
* Name:        add [static]
*
* Description: Addition of two polynomials modulo x^(2^log_dim) + 1.
*              Computation in a subring as used in SampleFz [GM18]
*              Works in normal and FFT representations
* 
* Arguments:   - fpr *out: polynomial to host the sum (initialized)
*              - const fpr *x: first polynomial factor
*              - const fpr *y: second polynomial factor
*              - const size_t log_dim: base-2 logarithm of the subring degree
**************************************************/
static inline void add(fpr *out, const fpr *x, const fpr *y, const size_t log_dim) {
  size_t n, i;

  n = (size_t)1 << log_dim;
  for (i = 0; i < n; i++) {
    out[i] = fpr_add(x[i], y[i]);
  }
}

/*************************************************
* Name:        sub [static]
*
* Description: Substraction of two polynomials modulo x^(2^log_dim) + 1. 
*              Computation in a subring as used in SampleFz [GM18].
*              Works in normal and FFT representations
* 
* Arguments:   - fpr *out: polynomial to host the difference (initialized)
*              - const fpr *x: first polynomial factor
*              - const fpr *y: second polynomial factor
*              - const size_t log_dim: base-2 logarithm of the subring degree
**************************************************/
static inline void sub(fpr *out, const fpr *x, const fpr *y, const size_t log_dim) {
  size_t n, i;

  n = (size_t)1 << log_dim;
  for (i = 0; i < n; i++) {
    out[i] = fpr_sub(x[i], y[i]);
  }
}

/*************************************************
* Name:        multiply [static]
*
* Description: Multiplication of two polynomials reduced
*              modulo x^(2^log_dim) + 1. Computation in 
*              a subring as used in SampleFz [GM18]
*              ONLY IN FFT REPRESENTATION
* 
* Arguments:   - fpr *out: polynomial to host the multiplication (initialized)
*              - const fpr *x: first polynomial factor
*              - const fpr *y: second polynomial factor
*              - const size_t log_dim: base-2 logarithm of the subring degree
**************************************************/
static inline void multiply(fpr *out, const fpr *x, const fpr *y, const size_t log_dim) {
  if (log_dim > 0) {
    size_t n, half_n, i;
    n = (size_t)1 << log_dim;

    // WARNING: Requires polynomial in FFT representation
    half_n = n >> 1;
    for (i = 0; i < half_n; i++) {
      fpr x_re, x_im, y_re, y_im;

      x_re = x[i];
      x_im = x[i + half_n];
      y_re = y[i];
      y_im = y[i + half_n];
      out[i] = fpr_sub(fpr_mul(x_re, y_re), fpr_mul(x_im, y_im));
      out[i + half_n] = fpr_add(fpr_mul(x_re, y_im), fpr_mul(x_im, y_re));
    }
  } else {
    out[0] = fpr_mul(x[0], y[0]);
  }
}

/*************************************************
* Name:        invert [static]
*
* Description: Invert a polynomial modulo x^(2^log_dim) + 1. 
*              Computation in a subring as used in SampleFz [GM18]
*              ONLY IN FFT REPRESENTATION
* 
* Arguments:   - fpr *out: polynomial to host the inverse (initialized)
*              - const fpr *x: polynomial to be inverted
*              - const size_t log_dim: base-2 logarithm of the subring degree
**************************************************/
static inline void invert(fpr *out, const fpr *x, const size_t log_dim) {
  if (log_dim > 0) {
    size_t n, half_n, i;
    n = (size_t)1 << log_dim;

    // WARNING: Requires polynomial in FFT representation
    half_n = n >> 1;
    for (i = 0; i < half_n; i++) {
      fpr x_re, x_im, tmp;

      x_re = x[i];
      x_im = x[i + half_n];
      tmp = fpr_add(fpr_mul(x_re, x_re), fpr_mul(x_im, x_im));

      out[i] = fpr_div(x_re, tmp);
      out[i + half_n] = fpr_neg(fpr_div(x_im, tmp));
    }
  } else {
    out[0] = fpr_inv(x[0]);
  }
}

/*************************************************
* Name:        conjugate [static]
*
* Description: Conjugate a polynomial modulo x^(2^log_dim) + 1. 
*              Computation in a subring as used in SampleFz [GM18]
*              ONLY IN FFT REPRESENTATION
* 
* Arguments:   - fpr *x: polynomial to be conjugated [in place] (initialized)
*              - size_t log_dim: base-2 logarithm of the subring degree
**************************************************/
static inline void conjugate(fpr *x, size_t log_dim) {
  size_t n, half_n, i;

  // WARNING: Requires polynomial in FFT representation
  n = (size_t)1 << log_dim;
  half_n = n >> 1;
  for (i = half_n; i < n; i++) {
    x[i] = fpr_neg(x[i]);
  }
}

/*************************************************
* Name:        split [static]
*
* Description: Split a polynomial of degree n into two polynomials
*              of degree n/2 by splitting the even and odd exponents. 
*              ONLY IN FFT REPRESENTATION
* 
* Arguments:   - fpr *even: polynomial to host the even exponents' coefficients (initialized)
*              - fpr *odd: polynomial to host the odd exponents' coefficients (initialized)
*              - const fpr *x: polynomial to be splitted
*              - size_t log_dim: base-2 logarithm of the subring degree
**************************************************/
static inline void split(fpr *even, fpr *odd, const fpr *x, const size_t log_dim) {
  size_t n, i;
  n = (size_t)1 << log_dim;

  // WARNING: Requires polynomial in FFT representation
  size_t half_n, quar_n;
  half_n = n >> 1;
  quar_n = half_n >> 1;

  even[0] = x[0];
  odd[0] = x[half_n];
  for (i = 0; i < quar_n; i++) {
    fpr a_re, a_im, b_re, b_im, t_re, t_im;
    a_re = x[(i << 1) + 0];
    a_im = x[(i << 1) + 0 + half_n];
    b_re = x[(i << 1) + 1];
    b_im = x[(i << 1) + 1 + half_n];

    t_re = fpr_add(a_re, b_re);
    t_im = fpr_add(a_im, b_im);
    even[i] = fpr_half(t_re);
    even[i + quar_n] = fpr_half(t_im);

    t_re = fpr_sub(a_re, b_re);
    t_im = fpr_sub(a_im, b_im);
    a_re = fpr_add(fpr_mul(t_re, fpr_gm_tab[((i + half_n) << 1) + 0]), fpr_mul(t_im, fpr_gm_tab[((i + half_n) << 1) + 1]));
    a_im = fpr_add(fpr_neg(fpr_mul(t_re, fpr_gm_tab[((i + half_n) << 1) + 1])), fpr_mul(t_im, fpr_gm_tab[((i + half_n) << 1) + 0]));
    odd[i] = fpr_half(a_re);
    odd[i + quar_n] = fpr_half(a_im);
  }
}

/*************************************************
* Name:        merge [static]
*
* Description: Merge two polynomials. Given even and odd
*              polynomials mod x^(2^(log_dim-1)) + 1, computes out = even(x^2) + x.odd(x^2)
*              mod x^(2^log_dim) + 1
*              ONLY IN FFT REPRESENTATION
* 
* Arguments:   - fpr *out: polynomial to host the merge (initialized)
*              - const fpr *even: polynomial with the even exponents' coefficients
*              - const fpr *odd: polynomial with the odd exponents' coefficients
*              - const size_t log_dim: base-2 logarithm of the subring degree
**************************************************/
static inline void merge(fpr *out, const fpr *even, const fpr *odd, const size_t log_dim) {
  size_t n, i;
  n = (size_t)1 << log_dim;

  // WARNING: Requires polynomial in FFT representation
  size_t half_n, quar_n;
  half_n = n >> 1;
  quar_n = half_n >> 1;

  out[0] = even[0];
  out[half_n] = odd[0];
  for (i = 0; i < quar_n; i++) {
    fpr a_re, a_im, b_re, b_im, t_re, t_im;
    a_re = even[i];
    a_im = even[i + quar_n];
    b_re = fpr_sub(fpr_mul(odd[i], fpr_gm_tab[((i + half_n) << 1) + 0]), fpr_mul(odd[i + quar_n], fpr_gm_tab[((i + half_n) << 1) + 1]));
    b_im = fpr_add(fpr_mul(odd[i], fpr_gm_tab[((i + half_n) << 1) + 1]), fpr_mul(odd[i + quar_n], fpr_gm_tab[((i + half_n) << 1) + 0]));

    t_re = fpr_add(a_re, b_re);
    t_im = fpr_add(a_im, b_im);
    out[(i << 1) + 0] = t_re;
    out[(i << 1) + 0 + half_n] = t_im;

    t_re = fpr_sub(a_re, b_re);
    t_im = fpr_sub(a_im, b_im);
    out[(i << 1) + 1] = t_re;
    out[(i << 1) + 1 + half_n] = t_im;
  }
}

/*************************************************
* Name:        _samplefz [static]
*
* Description: Sample a Gaussian element of Z[x]/(x^(2^log_dim) + 1)
*              of covariance real-valued polynomial f and center 
*              real-valued polynomial c, using the [GM18] SampleFz 
*              algorithm. Computation in a subring as used in [GM18]
*              Works in normal and FFT representations (using invert, multiply, conjugate, split, merge in adequate representation)
* 
* Arguments:   - fpr *res: polynomial to host the Gaussian sample (initialized)
*              - const fpr *f: covariance polynomial 
*              - const fpr *f: center polynomial 
*              - size_t log_dim: base-2 logarithm of the subring degree
**************************************************/
static void _samplefz(fpr *res, const fpr *f, const fpr *c, const size_t log_dim) {
  if (log_dim == 0) {
    res[0] = fpr_of(SampleZ(c[0], fpr_sqrt(f[0])));
  } else {
    size_t half_n = (size_t)1 << (log_dim - 1);
    fpr f_even[half_n], f_odd[half_n], c_even[half_n], c_odd[half_n], q_even[half_n], q_odd[half_n], tmp[half_n];

    // sampling begins
    split(f_even, f_odd, f, log_dim);
    split(c_even, c_odd, c, log_dim);

    _samplefz(q_odd, f_even, c_odd, log_dim-1);

    invert(tmp, f_even, log_dim-1);
    multiply(tmp, tmp, f_odd, log_dim-1);
    sub(c_odd, q_odd, c_odd, log_dim-1); // last usage of c_odd, this is now a temporary variable
    multiply(c_odd, tmp, c_odd, log_dim-1);
    add(c_even, c_even, c_odd, log_dim-1);

    conjugate(f_odd, log_dim-1);
    multiply(tmp, tmp, f_odd, log_dim-1);
    sub(tmp, f_even, tmp, log_dim-1); // last usage of f_even

    _samplefz(q_even, tmp, c_even, log_dim-1);

    // merge
    merge(res, q_even, q_odd, log_dim);
  }
}

/*************************************************
* Name:        poly_real_init
*
* Description: Initialize polynomial and set it to zero
*              This is strictly required before any operations 
*              are done with/on the polynomial.
* 
* Arguments:   - poly_real res: polynomial to be initialized
**************************************************/
void poly_real_init(poly_real res) {
  for (size_t i = 0; i < PARAM_N; i++) {
    res[i] = (fpr)PARAM_ZERO;
  }
}

/*************************************************
* Name:        poly_real_zero
*
* Description: Set polynomial to zero
* 
* Arguments:   - poly_real res: polynomial to be zeroized
**************************************************/
void poly_real_zero(poly_real res) {
  for (size_t i = 0; i < PARAM_N; i++) {
    res[i] = (fpr)PARAM_ZERO;
  }
}

/*************************************************
* Name:        poly_real_clear
*
* Description: Clears a polynomial and releases all associated memory. 
*              This is strictly required to avoid memory leaks and the 
*              polynomial must not be used again (unless reinitialized).
* 
* Arguments:   - poly_real res: polynomial to be cleared
**************************************************/
void poly_real_clear(poly_real res) {

}

/*************************************************
* Name:        poly_real_mul_scalar
*
* Description: Multiplication of a polynomial by a scalar
*              Works in normal and FFT representations
* 
* Arguments:   - poly_real res: polynomial to host the multiplication (initialized)
*              - const poly_real arg: first polynomial factor
*              - const coeff_real f: second scalar factor
**************************************************/
void poly_real_mul_scalar(poly_real res, const poly_real arg, const coeff_real f) {
  for (size_t i = 0; i < PARAM_N; i++) {
    res[i] = fpr_mul(arg[i], f);
  }
}

/*************************************************
* Name:        poly_real_add_constant
*
* Description: Add a scalar to the constant coefficient of a polynomial
*              ONLY IN FFT REPRESENTATION
* 
* Arguments:   - poly_real res: polynomial to be incremented (initialized)
*              - const coeff_real c: scalar summand
**************************************************/
void poly_real_add_constant(poly_real res, const coeff_real c) {
  // WARNING: Requires polynomial in FFT representation
  for (size_t i = 0; i < (PARAM_N >> 1); i++) {
    res[i] = fpr_add(res[i], c);
  }
}

/*************************************************
* Name:        poly_real_invert
*
* Description: Invert a real-valued polynomial modulo
*              x^PARAM_N + 1
* 
* Arguments:   - poly_real res: polynomial to host the inverse (initialized)
*              - const poly_real arg: polynomial to be inverted
**************************************************/
void poly_real_invert(poly_real res, const poly_real arg) {
  invert(res, arg, PARAM_LOG_N);
}

/*************************************************
* Name:        poly_real_mul
*
* Description: Multiplication of two polynomials reduced
*              modulo x^PARAM_N + 1
* 
* Arguments:   - poly_real res: polynomial to host the multiplication (initialized)
*              - const poly_real lhs: first polynomial factor
*              - const poly_real rhs: second polynomial factor
**************************************************/
void poly_real_mul(poly_real res, const poly_real lhs, const poly_real rhs) {
  multiply(res, lhs, rhs, PARAM_LOG_N);
}

/*************************************************
* Name:        poly_real_add
*
* Description: Add two polynomials 
* 
* Arguments:   - poly_real res: polynomial to host the sum (initialized)
*              - const poly_real lhs: first polynomial summand
*              - const poly_real rhs: second polynomial summand
**************************************************/
void poly_real_add(poly_real res, const poly_real lhs, const poly_real rhs) {
  add(res, lhs, rhs, PARAM_LOG_N);
}

/*************************************************
* Name:        poly_real_sub
*
* Description: Substract two polynomials 
* 
* Arguments:   - poly_real res: polynomial to host the difference (initialized)
*              - const poly_real lhs: first polynomial term
*              - const poly_real rhs: second polynomial term
**************************************************/
void poly_real_sub(poly_real res, const poly_real lhs, const poly_real rhs) {
  sub(res, lhs, rhs, PARAM_LOG_N);
}

/*************************************************
* Name:        poly_real_set_si
*
* Description: Set coefficient of x^n of a polynomial
*              Only used in normal representation (not FFT)
* 
* Arguments:   - poly_real arg: polynomial whose n-th coefficient is set (initialized)
*              - size_t n: degree of the coefficient to be set
*              - int64_t c: the new signed 64-bit int coefficient
**************************************************/
void poly_real_set_si(poly_real res, size_t n, int64_t c) {
  ASSERT_DEBUG(n < PARAM_N, "Illegal argument: cannot set coefficient of poly at given position.");
  res[n] = fpr_of(c);
}

/*************************************************
* Name:        poly_real_samplefz
*
* Description: Sample a Gaussian element of Z[x]/(x^PARAM_N + 1)
*              of covariance real-valued polynomial f and center 
*              real-valued polynomial c, using the [GM18] SampleFz 
*              algorithm
* 
* Arguments:   - poly_real res: polynomial to host the Gaussian sample (initialized)
*              - const poly_real f: covariance polynomial 
*              - const poly_real f: center polynomial
**************************************************/
void poly_real_samplefz(poly_real res, const poly_real f, const poly_real c) {
  _samplefz(res, f, c, PARAM_LOG_N);
}

/*************************************************
* Name:        poly_real_get_coeff_rounded
*
* Description: Get rounded coefficient of x^n of a polynomial
*              Only used in normal representation (not FFT)
* 
* Arguments:   - const poly_real arg: polynomial whose n-th coefficient is set (initialized)
*              - size_t n: degree of the coefficient to be set
*
* Returns the signed 64-bit integer rounded coefficient of x^n
**************************************************/
int64_t poly_real_get_coeff_rounded(const poly_real arg, size_t n) {
  ASSERT_DEBUG(n < PARAM_N, "Illegal argument: cannot get coefficient of poly at given position.");
  int64_t ret;
  ret = fpr_floor(fpr_add(arg[n], PARAM_ONE_HALF));
  return ret;
}
