#ifndef POLY_Q_SAMPLING_H
#define POLY_Q_SAMPLING_H

#include <stdint.h>
#include "params.h"
#include "arith.h"

void poly_q_uniform(poly_q pout, const uint8_t seed[SEED_BYTES], uint32_t domain_separator);
void poly_q_binomial(poly_q res, const uint8_t seed[SEED_BYTES], uint32_t cnt, uint32_t domain_separator);
void poly_q_vec_binomial(poly_q vec[PARAM_KH], const uint8_t seed[SEED_BYTES], uint32_t cnt, uint32_t domain_separator);

void poly_q_binary_fixed_weight(poly_q res, uint8_t state_in[STATE_BYTES]);

void poly_q_sample_gaussian_s2(poly_q res);
void poly_q_gaussian_coset_sL(poly_q res, const poly_q arg);
void poly_q_gaussian_coset_sH(poly_q res, const poly_q arg);

#endif /* POLY_Q_SAMPLING_H */
