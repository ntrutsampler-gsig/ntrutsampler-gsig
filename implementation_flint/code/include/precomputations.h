#ifndef PRECOMPUTATIONS_H
#define PRECOMPUTATIONS_H

#include <stdint.h>
#include "arith.h"

void fft_precomp_setup(void);
void fft_precomp_teardown(void);

void precompute_materials(poly_real g[PARAM_KH], poly_real f_prime[PARAM_KH], poly_real s_p, const poly_q f[PARAM_KH], const poly_q e[PARAM_KH], const poly_q e_cor_vec[PARAM_KH]);
int f_spectral_norm(const poly_q fi);
int e_sq_spectral_norm(poly_q e_cor_vec[PARAM_KH], const poly_q e[PARAM_KH]);

#endif /* PRECOMPUTATIONS_H */