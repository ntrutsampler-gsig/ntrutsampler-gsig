#include "sampling.h"
#include "fips202.h"
#include "arith.h"
#include "random.h"

/*************************************************
* Name:        perturbation_sampler
*
* Description: Gaussian perturbation sampler (PerturbationSampler, Algorithm 3.2)
* 
* Arguments:   - poly_q p1: polynomial to host top perturbation (initialized)
*              - poly_q *p2: array of polynomials to host bottom perturbation (initialized)
*              - const poly_real *g: array of polynomials, bottom sampling covariance
*              - const poly_real *f_prime: array of polynomials, center update vector
*              - const poly_real s_p: polynomial, Schur complement top sampling covariance
**************************************************/
void perturbation_sampler(
    poly_q p1, 
    poly_q p2[PARAM_KH], 
    const poly_real g[PARAM_KH],
    const poly_real f_prime[PARAM_KH],
    const poly_real s_p) {
  size_t i;
  poly_real c, tmp;

  // init polynomials
  poly_real_init(c);
  poly_real_init(tmp);

  poly_real_zero(c);
  for (i = 0; i < PARAM_KH; i++) {
    poly_real_zero(tmp);
    poly_q_samplefz(p2[i], g[i], tmp);
    poly_real_from_poly_q(tmp, p2[i]);
    poly_real_mul(tmp, tmp, f_prime[i]);
    poly_real_add(c, c, tmp);
  }
  poly_q_samplefz(p1, s_p, c);

  //  clean up polynomials
  poly_real_clear(c);
  poly_real_clear(tmp);
}

/*************************************************
* Name:        gadget_sampler
*
* Description: Gadget sampling of zL and zH such that [1 | q_L.g_{H,t}^T][z_L | z_H] = w
*              (Steps 3 to 10 of Sampler)
* 
* Arguments:   - poly_q zL: polynomial to host z_L (initialized)
*              - poly_q *zH: array of polynomials to host z_H (initialized)
*              - const poly_q w: polynomial corrected syndrome w
*              - const poly_q tag: polynomial, tag
*              - const poly_q taginv: polynomial, tag^{-1} mod bH
*              - const poly_q *f: array of polynomials containing the f_i 
*              - const poly_q *finv_bH: array of polynomials containing the f_i^{-1} mod bH 
**************************************************/
void gadget_sampler(
    poly_q zL,
    poly_q zH[PARAM_KH],
    const poly_q w,
    const poly_q tag,
    const poly_q taginv,
    const poly_q f[PARAM_KH],
    const poly_q finv_bH[PARAM_KH]) {
  size_t i;
  poly_q tmp_poly, tmp_center;

  // init polynomials
  poly_q_init(tmp_poly);
  poly_q_init(tmp_center);

  // sampling z_L
  poly_q_set(tmp_poly, w);
  poly_q_mod_qL(tmp_center, tmp_poly);
  poly_q_gaussian_coset_sL(zL, tmp_center);
  poly_q_sub(tmp_center, w, zL); // tmp_center used as tmp variable
  poly_q_div_qL(tmp_poly, tmp_center);

  // sampling zH
  for (i = 0; i < PARAM_KH; i++) {
    poly_q_mul(tmp_center, taginv, tmp_poly);
    poly_q_mul(tmp_center, tmp_center, finv_bH[i]);
    poly_q_mod_bH(tmp_center, tmp_center);
    poly_q_gaussian_coset_sH(zH[i], tmp_center);
    poly_q_mul(tmp_center, tag, zH[i]); // tmp_center used as tmp variable
    poly_q_mul(tmp_center, tmp_center, f[i]); 
    poly_q_sub(tmp_center, tmp_poly, tmp_center);
    poly_q_div_bH(tmp_poly, tmp_center);
  }

  // clean up polynomials
  poly_q_clear(tmp_poly);
  poly_q_clear(tmp_center);
}

/*************************************************
* Name:        ntru_tsampler
*
* Description: NTRU-TSampler
* 
* Arguments:   - poly_q_vec_d *v1: array of polynomial vectors to host top preimage (initialized)
*              - poly_q_vec_d *v2: array of polynomial vectors to host bottom preimage (initialized)
*              - const poly_q_mat_d_d *R: array of polynomial matrices, secret key R 
*              - const poly_q_mat_d_d A: polynomial matrix, public A'
*              - const poly_q_mat_d_d *B: array of polynomial matrices, public B
*              - const poly_q_vec_d u: polynomial vector, public u
*              - const poly_q tag: polynomial, tag
*              - const poly_q taginv: polynomial, tag^{-1} mod bH
*              - const poly_real_mat_2d_2d *S: array of polynomial matrices, Schur complements for sampling
**************************************************/
void ntru_tsampler(
    poly_q v1, 
    poly_q v2[PARAM_KH], 
    const poly_q h[PARAM_KH], 
    const poly_q u, 
    const poly_q tag, 
    const poly_q taginv,
    const isk_t *isk) {
  size_t i;
  poly_q p1;
  poly_q p2[PARAM_KH];
  poly_q w, zL, zH[PARAM_KH];
  int64_t bexpi;

  // init polynomials
  poly_q_init(p1);
  poly_q_init(w);
  poly_q_init(zL);
  for (i = 0; i < PARAM_KH; i++) {
    poly_q_init(p2[i]);
    poly_q_init(zH[i]);
  }

  // perturbation sampler
  perturbation_sampler(p1, p2, isk->g, isk->f_prime, isk->s_p);

  // w = u - p1 + (h - qL.tag.gH)^T.p2
  poly_q_sub(w, u, p1);
  bexpi = PARAM_QL;
  for (i = 0; i < PARAM_KH; i++) {
    poly_q_mul_scalar(zL, tag, bexpi); // using zL as temporary variable
    poly_q_sub(zL, h[i], zL);
    poly_q_mul(zL, zL, p2[i]);
    poly_q_add(w, w, zL);
    bexpi *= PARAM_BH;
  }

  // Gadget sampling
  gadget_sampler(zL, zH, w, tag, taginv, isk->f, isk->finv_bH);

  // v1 = p1 + zL + sum_j e[j].zH[j] and v2[i] = p2[i] + f[i].zH[i]
  poly_q_add(v1, p1, zL);
  for (i = 0; i < PARAM_KH; i++) {
    poly_q_mul(zL, isk->e[i], zH[i]); // using zL as temporary variable
    poly_q_add(v1, v1, zL);
    poly_q_mul(zL, isk->f[i], zH[i]); // using zL as temporary variable
    poly_q_add(v2[i], p2[i], zL);
  }

  // clean up vectors
  poly_q_clear(p1);
  poly_q_clear(w);
  poly_q_clear(zL);
  for (i = 0; i < PARAM_KH; i++) {
    poly_q_clear(p2[i]);
    poly_q_clear(zH[i]);
  }
}