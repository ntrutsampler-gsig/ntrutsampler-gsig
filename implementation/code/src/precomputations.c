#include "precomputations.h"
#include "arith.h"
#include "fpr.h"
#include "fft.h"

/*************************************************
* Name:        f_sq_spectral_norm
*
* Description: Computes the spectral norm of polynomial fi
*              as ||sigma(fi)||_infty
*
* Arguments:   - const poly_q fi: polynomial hosting fi
* 
* Returns an int whether or not the norm exceeds PARAM_BF
**************************************************/
int f_sq_spectral_norm(const poly_q fi) {
  size_t j;
  int is_invalid = 0;
  fpr sq_abs_val, vec[PARAM_N];

  // Compute FFT (for sigma embedding)
  for (j = 0; j < PARAM_N; j++) {
    vec[j] = poly_q_get_coeff_centered(fi, j);
  }
  FFT(vec, PARAM_LOG_N);

  // Check infinity norm
  for (j = 0; j < PARAM_N / 2; j++) {
    sq_abs_val = fpr_add(fpr_mul(vec[j], vec[j]), fpr_mul(vec[j+PARAM_N/2], vec[j+PARAM_N/2]));
    if (fpr_lt(PARAM_BFSQ, sq_abs_val)) {
      is_invalid = 1;
      goto f_spectral_norm_cleanup;
    }
  }

f_spectral_norm_cleanup:
  return is_invalid;
}

/*************************************************
* Name:        e_sq_spectral_norm
*
* Description: Computes the spectral norm of polynomial e*.e
*              as ||sigma(e*.e)||_infty. Stores the e_i*.e_i in 
*              e_cor_vec
*
* Arguments:   - poly_q *e_cor_vec: array of polynomial hosting e_i*.e_i (initialized)
*              - const poly_q *e: array of polynomial hosting e
* 
* Returns an int whether or not the norm exceeds PARAM_BESQ
**************************************************/
int e_sq_spectral_norm(poly_q e_cor_vec[PARAM_KH], const poly_q e[PARAM_KH]) {
  size_t j;
  int is_invalid = 0;
  poly_q corr;
  fpr sq_abs_val, vec[PARAM_N];

  // init flint objects 
  poly_q_init(corr);

  // Computing correlation e*.e (modulo q, but small elements so there is no wraparound)
  poly_q_zero(corr);
  for (j = 0; j < PARAM_KH; j++) {
    poly_q_conjugate(e_cor_vec[j], e[j]);
    poly_q_mul(e_cor_vec[j], e_cor_vec[j], e[j]);
    poly_q_add(corr, corr, e_cor_vec[j]);
  }

  // Compute FFT (for sigma embedding)
  for (j = 0; j < PARAM_N; j++) {
    vec[j] = poly_q_get_coeff_centered(corr, j);
  }
  FFT(vec, PARAM_LOG_N);

  // Check infinity norm
  for (j = 0; j < PARAM_N / 2; j++) {
    sq_abs_val = fpr_add(fpr_mul(vec[j], vec[j]), fpr_mul(vec[j+PARAM_N/2], vec[j+PARAM_N/2]));
    if (fpr_lt(PARAM_BE_POW_4, sq_abs_val)) {
      is_invalid = 1;
      goto e_sq_spectral_norm_cleanup;
    }
  }

e_sq_spectral_norm_cleanup:
  // cleanup flint objects
  poly_q_clear(corr);

  return is_invalid;
}

/*************************************************
* Name:        precompute_materials
*
* Description: Computes Schur complements of covariance for perturbation 
*              sampling and center update from the issuer's secret key
*
* Arguments:   - poly_real *g: 
**************************************************/
void precompute_materials(
    poly_real g[PARAM_KH], 
    poly_real f_prime[PARAM_KH], 
    poly_real s_p, 
    const poly_q f[PARAM_KH], 
    const poly_q e[PARAM_KH], 
    const poly_q e_cor_vec[PARAM_KH]) {
  
  size_t i;
  poly_q tmp;
  poly_real tmp1, tmp2;

  // init polynomials
  poly_q_init(tmp);
  poly_real_init(tmp1);
  poly_real_init(tmp2);

  poly_real_zero(s_p);
  for (i = 0; i < PARAM_KH; i++) {
    // Computing g[i] = s_2² - s_H².f_i*.f_i
    poly_q_conjugate(tmp, f[i]);
    poly_real_from_poly_q(g[i], f[i]);
    poly_real_from_poly_q(tmp1, tmp);
    poly_real_mul(g[i], tmp1, g[i]);
    poly_real_mul_scalar(g[i], g[i], PARAM_NEG_SHSQ);
    poly_real_add_constant(g[i], PARAM_S2SQ);

    // Computing -s_H².g[i]^{-1} in f_prime[i]
    poly_real_invert(f_prime[i], g[i]);
    poly_real_mul_scalar(f_prime[i], f_prime[i], PARAM_NEG_SHSQ);

    // Converting e_i*.e_i to poly_real for s_p
    poly_real_from_poly_q(tmp2, e_cor_vec[i]);
    poly_real_mul(tmp2, tmp2, f_prime[i]); // contains -s_H².e_i*.e_i.(s_2² - s_H².f_i*.f_i)^{-1}
    poly_real_add(s_p, s_p, tmp2);

    // Finishing computation of f_prime[i] (tmp1 holds f_i*)
    poly_real_from_poly_q(tmp2, e[i]);
    poly_real_mul(tmp2, tmp1, tmp2); // contains e_i.f_i*
    poly_real_mul(f_prime[i], f_prime[i], tmp2);
  }

  // Finishing computation of s_p
  poly_real_mul_scalar(s_p, s_p, PARAM_S2SQ);
  poly_real_add_constant(s_p, PARAM_S1SQ_SLSQ);

  // clean up polynomials
  poly_q_clear(tmp);
  poly_real_clear(tmp1);
  poly_real_clear(tmp2);
}