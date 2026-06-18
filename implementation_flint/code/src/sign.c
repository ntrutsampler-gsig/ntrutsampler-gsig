#include "arith.h"
#include "randombytes.h"
#include "sampling.h"
#include "sign.h"
#include "poly_q_sampling.h"
#include "fips202.h"
#include "precomputations.h"

/*************************************************
* Name:        issuer_keys_init
*
* Description: Initialize keys of sampler
*              by calling Flint initialization
*
* Arguments:   - ipk_t *ipk: pointer to public key structure
*              - isk_t *isk: pointer to secret key structure
**************************************************/
void issuer_keys_init(ipk_t *ipk, isk_t *isk) {
  size_t i;
  for (i = 0;  i < PARAM_KH; i++) {
    poly_q_init(ipk->h[i]);

    poly_q_init(isk->e[i]);
    poly_q_init(isk->f[i]);
    poly_q_init(isk->finv_bH[i]);
    poly_real_init(isk->g[i]);
    poly_real_init(isk->f_prime[i]);
  }
  poly_real_init(isk->s_p);
}

/*************************************************
* Name:        issuer_keys_clear
*
* Description: Clear keys of sampler
*              by calling Flint clean up
*
* Arguments:   - ipk_t *ipk: pointer to public key structure
*              - isk_t *isk: pointer to secret key structure
**************************************************/
void issuer_keys_clear(ipk_t *ipk, isk_t *isk) {
  size_t i;
  for (i = 0;  i < PARAM_KH; i++) {
    poly_q_clear(ipk->h[i]);

    poly_q_clear(isk->e[i]);
    poly_q_clear(isk->f[i]);
    poly_q_clear(isk->finv_bH[i]);
    poly_real_clear(isk->g[i]);
    poly_real_clear(isk->f_prime[i]);
  }
  poly_real_clear(isk->s_p);
}

/*************************************************
* Name:        usk_init
*
* Description: Initialize structure to host the user secret key
*              by calling Flint initialization
*
* Arguments:   - usk_t *usk: pointer to user key structure
**************************************************/
void usk_init(usk_t *usk) {
  poly_q_init(usk->tag);
  for (size_t i = 0; i < PARAM_KH; i++) {
    poly_q_init(usk->v2[i]);
  }
  poly_q_init(usk->v3);
}

/*************************************************
* Name:        usk_clear
*
* Description: Clear structure that hosts the user secret key
*              by calling Flint clean up
*
* Arguments:   - usk_t *usk: pointer to user key structure
**************************************************/
void usk_clear(usk_t *usk) {
  poly_q_clear(usk->tag);
  for (size_t i = 0; i < PARAM_KH; i++) {
    poly_q_clear(usk->v2[i]);
  }
  poly_q_clear(usk->v3);
}

/*************************************************
* Name:        issuer_keygen
*
* Description: Generates public and private key of the issuer
*
* Arguments:   - ipk_t *ipk: pointer to issuer public key structure (initialized)
*              - isk_t *isk: pointer to issuer secret key structure (initialized)
**************************************************/
void issuer_keygen(ipk_t *ipk, isk_t *isk) {
  uint8_t root_seed[SEED_BYTES], seeds[SEED_BYTES*2];
  uint8_t *public_seed = seeds, *secret_seed = &seeds[SEED_BYTES];
  poly_q e_correlation_vec[PARAM_KH];
  uint32_t kappa;
  size_t i;

  // generate random seed(s)
  randombytes(root_seed, SEED_BYTES);
  sha3_512(seeds, root_seed, SEED_BYTES);
#if SEED_BYTES != 32
#error "SEED_BYTES must be 32."
#endif

  // init vector
  for (i = 0; i < PARAM_KH; i++) {
    poly_q_init(e_correlation_vec[i]);
  }
  
  // sample f from B_1 (binomial)
  kappa = 0;
  for (i = 0; i < PARAM_KH; i++) {
    do {
      poly_q_binomial(isk->f[i], secret_seed, kappa++, DOMAIN_SEPARATOR_F);
    } while(f_spectral_norm(isk->f[i]));
  }

  // sample e from B_1 (binomial) and computes 
  kappa = 0;
  do {
    poly_q_vec_binomial(isk->e, secret_seed, kappa++, DOMAIN_SEPARATOR_E);
  } while(e_sq_spectral_norm(e_correlation_vec, isk->e));

  // compute g = [g_i]_i = [s_2² - s_H².f_i.f_i*]_i
  //         f'= [-s_H².e_i.f_i*.g_i^{-1}]_i
  //         s_p = s_1² - s_L² - s_H².s_2².sum_i e_i.e_i*.g_i^{-1}
  precompute_materials(isk->g, isk->f_prime, isk->s_p, isk->f, isk->e, e_correlation_vec);

  // compute inverses of f_i modulo q first, compute h_i = e_i / f_i mod q, and reduce mod b_H to get inverse mod b_H
  // Using e_correlation_vec[0] as temporary variable
  for (i = 0; i < PARAM_KH; i++) {
    poly_q_invert_mod_q_bH(e_correlation_vec[0], isk->finv_bH[i], isk->f[i]);
    poly_q_mul(ipk->h[i], isk->e[i], e_correlation_vec[0]);
    //poly_q_mod_bH(isk->finv_bH[i], e_correlation_vec[0]);
  }

  // append public seed to pk for extending u, a3
  for (i = 0; i < SEED_BYTES; i++) {
    ipk->seed[i] = public_seed[i];
  }

  // clean up
  for (i = 0; i < PARAM_KH; i++) {
    poly_q_clear(e_correlation_vec[i]);
  }
}

/*************************************************
* Name:        tag_gen
*
* Description: Computes the tag from the state
*
* Arguments:   - poly_q tag: polynomial hosting the tag used for commitment/signature (initialized)
*              - uint8_t *state: pointer to signer's state byte array (allocated STATE_BYTES bytes)
**************************************************/
void tag_gen(poly_q tag, uint8_t state[STATE_BYTES]) {
  size_t i;

  // compute tag from state
  poly_q_binary_fixed_weight(tag, state);

  // increment state
#if (STATE_BYTES%4) != 0
#error "STATE_BYTES must be multiple of 4"
#endif
  uint64_t stateinc = ((uint64_t)*(uint32_t*)state) + 1; // cast byte array to uint32, then to uint64, then increment
  *(uint32_t*)state = (uint32_t) stateinc;
  uint64_t carry = stateinc >> 32;
  for (i = 1; i < STATE_BYTES/4; i++)
  {
    stateinc = ((uint64_t)*(uint32_t*)&state[i*4]) + carry; // cast byte array to uint32, then to uint64, then add carry
    *(uint32_t*)&state[i*4] = (uint32_t) stateinc;
    carry = stateinc >> 32;
  }
}

/*************************************************
* Name:        issuer_sign
*
* Description: Computes the user secret key 
*
* Arguments:   - usk_t *usk: pointer to user key structure (initialized)
*              - uint8_t *state: pointer to signer's state byte array (allocated STATE_BYTES bytes)
*              - const isk_t *isk: pointer to issuer secret key structure
*              - const ipk_t *ipk: pointer to issuer public key structure
**************************************************/
void issuer_sign(usk_t *usk, uint8_t state[STATE_BYTES], const isk_t *isk, const ipk_t *ipk) {
  size_t i;
  poly_q u, a3, tmp, taginv, v1;
  uint64_t norm2sq_v1, norm2sq_v2_v3;

  // init polynomials
  poly_q_init(u);
  poly_q_init(a3);
  poly_q_init(tmp);
  poly_q_init(taginv);
  poly_q_init(v1);

  // expand uniform elements a3, u
  poly_q_uniform(a3, ipk->seed, DOMAIN_SEPARATOR_A3);
  poly_q_uniform(u, ipk->seed, DOMAIN_SEPARATOR_U);

  // compute tag and update state
  tag_gen(usk->tag, state);

  // invert tag modulo PARAM_BH in m (used as tmp variable)
  poly_q_invert_mod_bH(taginv, usk->tag);

reject_preimage:
  // sample v3 from discrete Gaussian
  poly_q_sample_gaussian_s2(usk->v3); // probabilistic

  // compute u - a3.v3
  poly_q_mul(tmp, a3, usk->v3);
  poly_q_sub(tmp, u, tmp);

  // call NTRU-TSampler, output in v1, usk->v2
  ntru_tsampler(v1, usk->v2, ipk->h, tmp, usk->tag, taginv, isk); // probabilistic (tmp holds u - a3.v3)

  // check l2 norms before outputting preimage
  norm2sq_v1 = poly_q_sq_norm2(v1);
  norm2sq_v2_v3 = poly_q_sq_norm2(usk->v3);
  for (i = 0; i < PARAM_KH; i++) {
    norm2sq_v2_v3 += poly_q_sq_norm2(usk->v2[i]);
  }
  
  if((norm2sq_v1 > PARAM_B1SQ) || (norm2sq_v2_v3 > PARAM_B2SQ)) {
    goto reject_preimage;
  }

  // clean up polynomials
  poly_q_clear(u);
  poly_q_clear(a3);
  poly_q_clear(tmp);
  poly_q_clear(taginv);
  poly_q_clear(v1);
}

/*************************************************
* Name:        user_verify
*
* Description: Verifies preimage is correct (before ZKP in group signature)
*
* Arguments:   - const usk_t *usk: pointer to the user key structure 
*              - const ipk_t *ipk: pointer to issuer public key structure
* 
* Returns 1 if signature could be verified correctly and 0 otherwise
**************************************************/
int user_verify(const usk_t *usk, const ipk_t *ipk) {
  size_t i;
  int64_t bexpi;
  poly_q u, a3, v1;
  uint64_t norm2sq_v1, norm2sq_v2_v3;
  int64_t tag_weight;

  // init polynomials
  poly_q_init(u);
  poly_q_init(a3);
  poly_q_init(v1);

  // expand uniform elements a3, u
  poly_q_uniform(a3, ipk->seed, DOMAIN_SEPARATOR_A3);
  poly_q_uniform(u, ipk->seed, DOMAIN_SEPARATOR_U);

  // compute v1 = u + (h - qL.t.gH)^T.v2 - a3.v3
  poly_q_mul(v1, a3, usk->v3); // using v1 as tmp variable
  poly_q_sub(v1, u, v1); 
  bexpi = PARAM_QL;
  for (i = 0; i < PARAM_KH; i++) {
    poly_q_mul_scalar(a3, usk->tag, bexpi); // using a3 as tmp variable
    poly_q_sub(a3, ipk->h[i], a3);
    poly_q_mul(a3, a3, usk->v2[i]);
    poly_q_add(v1, v1, a3);
    bexpi *= PARAM_BH;
  }

  // check l2 norms before outputting preimage
  norm2sq_v1 = poly_q_sq_norm2(v1);
  norm2sq_v2_v3 = poly_q_sq_norm2(usk->v3);
  for (i = 0; i < PARAM_KH; i++) {
    norm2sq_v2_v3 += poly_q_sq_norm2(usk->v2[i]);
  }
  tag_weight = poly_q_weight(usk->tag); // returns -1 if polynomial is not in I(x^{k_pi}), else the number of ones
  // clean up polynomials
  poly_q_clear(u);
  poly_q_clear(a3);
  poly_q_clear(v1);

  // return check
  return (norm2sq_v1 <= PARAM_B1SQ) && (norm2sq_v2_v3 <= PARAM_B2SQ) && (tag_weight == PARAM_W);
}