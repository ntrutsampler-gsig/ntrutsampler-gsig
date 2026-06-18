#ifndef SAMPLING_H
#define SAMPLING_H

#include <stdint.h>
#include "params.h"
#include "arith.h"
#include "sign.h"

void perturbation_sampler(
    poly_q p1, 
    poly_q p2[PARAM_KH], 
    const poly_real g[PARAM_KH],
    const poly_real f_prime[PARAM_KH],
    const poly_real s_p);
void gadget_sampler(
    poly_q zL,
    poly_q zH[PARAM_KH],
    const poly_q w,
    const poly_q tag,
    const poly_q taginv,
    const poly_q f[PARAM_KH],
    const poly_q finv_bH[PARAM_KH]);
void ntru_tsampler(
    poly_q v1, 
    poly_q v2[PARAM_KH], 
    const poly_q h[PARAM_KH], 
    const poly_q u, 
    const poly_q tag, 
    const poly_q taginv,
    const isk_t *isk);

#endif /* SAMPLING_H */
