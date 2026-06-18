#include "poly_q_sampling.h"
#include "fips202.h"
#include "arith.h"
#include "random.h"

/*************************************************
* Name:        poly_q_uniform
*
* Description: Sample a uniformly random polynomial modulo
*              PARAM_Q deterministically from a seed.
* 
* Arguments:   - poly_q pout: output uniform polynomial (initialized)
*              - const uint8_t *seed: pointer to byte array containing seed (allocated SEED_BYTES bytes)
*              - uint32_t domain_separator: domain separator for XOF
**************************************************/
void poly_q_uniform(poly_q pout, const uint8_t seed[SEED_BYTES], uint32_t domain_separator) {
  uint8_t output[SHAKE128_RATE * 2];
  keccak_state state;
  size_t k,cnt,off,bytecnt;
  shake128_init(&state);
  shake128_absorb(&state, seed, SEED_BYTES);
  shake128_absorb(&state, (const uint8_t*)&domain_separator, sizeof(uint32_t));
  shake128_finalize(&state);
  shake128_squeezeblocks(output, 2, &state);
  bytecnt = 2*SHAKE128_RATE;

  cnt = 0;
  off = 0;
  while (cnt < PARAM_N) {
#if PARAM_Q_BITLEN > 18
#error "PARAM_Q_BITLEN too big for uniform sampling."
#else
    // idea: take 7 bytes, partition into 3 chunks of 17 bits, potentially ignore the MSBs, perform rejection sampling
    if (bytecnt < 7) {
      for (k = 0; k < bytecnt; k++) {
        output[k] = output[off++];
      }
      shake128_squeezeblocks(&output[bytecnt], 1, &state);
      off = 0;
      bytecnt += SHAKE128_RATE;
    }
    uint64_t tmp7byte = output[off] | ((uint64_t)output[off+1] << 8) | ((uint64_t)output[off+2] << 16) | ((uint64_t)output[off+3] << 24) | ((uint64_t)output[off+4] << 32) | ((uint64_t)output[off+5] << 40) | ((uint64_t)output[off+6] << 48);
    uint64_t tmp = tmp7byte & ((1<<PARAM_Q_BITLEN)-1);
    if (tmp < (coeff_q)PARAM_Q) {
      poly_q_set_coeff(pout, cnt++, tmp);
      if (cnt == PARAM_N) {
        break;
      }
    }
    tmp7byte = tmp7byte >> PARAM_Q_BITLEN;
    tmp = tmp7byte & ((1<<PARAM_Q_BITLEN)-1);
    if (tmp < (coeff_q)PARAM_Q) {
      poly_q_set_coeff(pout, cnt++, tmp);
      if (cnt == PARAM_N) {
        break;
      }
    }
    tmp = (tmp7byte >> PARAM_Q_BITLEN) & ((1<<PARAM_Q_BITLEN)-1);
    if (tmp < (coeff_q)PARAM_Q) {
      poly_q_set_coeff(pout, cnt++, tmp);
    }
  
    off += 7;
    bytecnt -= 7;
#if PARAM_Q_BITLEN < 17
#warning "PARAM_Q_BITLEN maybe unsuitable for efficient uniform sampling."
#endif
#endif
  }
}

/*************************************************
* Name:        poly_q_binomial
*
* Description: Sample a centered binomial polynomial with binomial 
*              parameter 1 deterministically from a seed.
* 
* Arguments:   - poly_q res: output binomial polynomial (initialized)
*              - const uint8_t *seed: pointer to byte array containing seed (allocated SEED_BYTES bytes)
*              - uint32_t cnt: repetition domain separator for XOF
*              - uint32_t domain_separator: domain separator for XOF
**************************************************/
void poly_q_binomial(poly_q res, const uint8_t seed[SEED_BYTES], uint32_t cnt, uint32_t domain_separator) {
#if (PARAM_N%64) != 0
#error "PARAM_N must be divisible by 64"
#endif
  uint64_t output[PARAM_N*2/64]; // 2 bits per coefficient
  uint64_t coef_lsb[PARAM_N/64];
  uint64_t coef_sign[PARAM_N/64];
  keccak_state state;
  size_t i;

  shake256_init(&state);
  shake256_absorb(&state, seed, SEED_BYTES);
  shake256_absorb(&state, (const uint8_t*)&cnt, sizeof(uint32_t));
  shake256_absorb(&state, (const uint8_t*)&domain_separator, sizeof(uint32_t));
  shake256_finalize(&state);
  shake256_squeeze((uint8_t*)output, PARAM_N*2/8, &state);

  for (i = 0; i < PARAM_N/64; i++) {
    coef_lsb[i] = output[2*i] ^ output[2*i+1];
    coef_sign[i] = output[2*i] & output[2*i+1];
  }
  for (i = 0; i < PARAM_N; i++) {
    poly_q_set_coeff(res, i, (int32_t)((coef_lsb[i/64] >> (i%64))&1) + (int32_t)(((coef_sign[i/64] >> ((i%64))) << 1)&2) - 1);
    // we have for sign||lsb either 00 (->-1) or 01 (->0) or 10 (->1), so we reconstruct this and subtract one
  }
}

/*************************************************
* Name:        poly_q_vec_binomial
*
* Description: Sample a centered binomial polynomial vector of
*              size PARAM_KH with binomial parameter 1 
*              deterministically from a seed.
* 
* Arguments:   - poly_q *vec: pointer output binomial polynomial vector (initialized)
*              - const uint8_t *seed: pointer to byte array containing seed (allocated SEED_BYTES bytes)
*              - uint32_t cnt: repetition domain separator for XOF
*              - uint32_t domain_separator: domain separator for XOF
**************************************************/
void poly_q_vec_binomial(poly_q vec[PARAM_KH], const uint8_t seed[SEED_BYTES], uint32_t cnt, uint32_t domain_separator) {
#if (PARAM_N%64) != 0
#error "PARAM_N must be divisible by 64"
#endif
  uint64_t output[PARAM_N*2/64]; // 2 bits per coefficient
  uint64_t coef_lsb[PARAM_N/64];
  uint64_t coef_sign[PARAM_N/64];
  keccak_state state;
  size_t i,j;

  for (i = 0; i < PARAM_KH; i++) {
    shake256_init(&state);
    shake256_absorb(&state, seed, SEED_BYTES);
    shake256_absorb(&state, (const uint8_t*)&cnt, sizeof(uint32_t));
    shake256_absorb(&state, (const uint8_t*)&domain_separator, sizeof(uint32_t));
    shake256_absorb(&state, (const uint8_t*) &i, sizeof(i));
    shake256_finalize(&state);
    shake256_squeeze((uint8_t*)output, PARAM_N*2/8, &state);

    for (j = 0; j < PARAM_N/64; j++) {
      coef_lsb[j] = output[2*j] ^ output[2*j+1];
      coef_sign[j] = output[2*j] & output[2*j+1];
    }
    for (j = 0; j < PARAM_N; j++) {
      poly_q_set_coeff(vec[i], j, (int32_t)((coef_lsb[j/64] >> (j%64))&1) + (int32_t)(((coef_sign[j/64] >> ((j%64))) << 1)&2) - 1);
      // we have for sign||lsb either 00 (->-1) or 01 (->0) or 10 (->1), so we reconstruct this and subtract one
    }
  }
}

/*************************************************
* Name:        poly_q_binary_fixed_weight
*
* Description: Sample binary polynomial id with fixed Hamming weight PARAM_W
*              in R_pi and generates res = id(X^k_pi)
* 
* Arguments:   - poly_q res: output uniform binary polynomial with fixed Hamming weight (initialized)
*              - uint8_t *state_in: state from which to expand the polynomial (allocated STATE_BYTES bytes)
**************************************************/
void poly_q_binary_fixed_weight(poly_q res, uint8_t state_in[STATE_BYTES]) {
  unsigned int i, b, pos = 0;
  uint8_t buf[SHAKE256_RATE];
  keccak_state state;

  shake256_absorb_once(&state, state_in, STATE_BYTES);
  shake256_squeezeblocks(buf, 1, &state);

  for (i = 0; i < PARAM_N; i++) {
    poly_q_set_coeff(res, i, 0);
  }
  for (i = PARAM_N_PI - PARAM_W; i < PARAM_N_PI; i++) {
    do {
      if (pos >= SHAKE256_RATE) {
        shake256_squeezeblocks(buf, 1, &state);
        pos = 0;
      }
      b = buf[pos++];
    } while (b > i);

    poly_q_set_coeff(res, i * PARAM_K_PI, poly_q_get_coeff(res, b * PARAM_K_PI));
    poly_q_set_coeff(res, b * PARAM_K_PI, 1);
  }
}
