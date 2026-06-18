#ifndef SIGN_H
#define SIGN_H

#include <stdint.h>
#include "params.h"
#include "arith.h"

typedef struct {
  poly_q e[PARAM_KH];
  poly_q f[PARAM_KH];

  poly_q finv_bH[PARAM_KH];
  poly_real g[PARAM_KH];
  poly_real f_prime[PARAM_KH];
  poly_real s_p;
} isk_t;

typedef struct {
  poly_q h[PARAM_KH];
  uint8_t seed[SEED_BYTES];
} ipk_t;

typedef struct {
  poly_q tag;
  poly_q v2[PARAM_KH];
  poly_q v3;
} usk_t;

void issuer_keys_init(ipk_t *ipk, isk_t *isk);
void issuer_keys_clear(ipk_t *ipk, isk_t *isk);
void usk_init(usk_t *usk);
void usk_clear(usk_t *usk);

void issuer_keygen(ipk_t *ipk, isk_t *isk);
void tag_gen(poly_q tag, uint8_t state[STATE_BYTES]);
void issuer_sign(usk_t *usk, uint8_t state[STATE_BYTES], const isk_t *isk, const ipk_t *ipk);
int user_verify(const usk_t *usk, const ipk_t *ipk);

#endif /* SIGN_H */
