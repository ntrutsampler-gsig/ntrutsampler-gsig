#include "bench_sign.h"

#include "arith.h"
#include "sign.h"
#include "randombytes.h"
#include "random.h"

double issuer_keygen_bench(timer* t) {
  double time;
  isk_t isk;
  ipk_t ipk;
  /* ----------- BEGIN: Code under measurement --------- */
  start_timer(t);
  issuer_keys_init(&ipk, &isk);
  issuer_keygen(&ipk, &isk);
  time = stop_timer(t);
  /* ----------- END: Code under measurement ----------- */
  issuer_keys_clear(&ipk, &isk);
  return time;
}

double issuer_sign_bench(timer* t) {
  double time;
  isk_t isk;
  ipk_t ipk;
  usk_t usk;
  uint8_t state[STATE_BYTES];

  issuer_keys_init(&ipk, &isk);
  usk_init(&usk);
  issuer_keygen(&ipk, &isk);
  randombytes(state, STATE_BYTES);
  /* ----------- BEGIN: Code under measurement --------- */
  start_timer(t);
  issuer_sign(&usk, state, &isk, &ipk);
  time = stop_timer(t);
  /* ----------- END: Code under measurement ----------- */
  if (!user_verify(&usk, &ipk)) {
    printf("FATAL ERROR: benchmarked user key is not valid\n");
  }
  usk_clear(&usk);
  issuer_keys_clear(&ipk, &isk);
  return time;
}

double user_verify_valid_bench(timer* t) {
  double time;
  isk_t isk;
  ipk_t ipk;
  usk_t usk;
  uint8_t state[STATE_BYTES];

  issuer_keys_init(&ipk, &isk);
  usk_init(&usk);
  issuer_keygen(&ipk, &isk);
  randombytes(state, STATE_BYTES);
  issuer_sign(&usk, state, &isk, &ipk);
  /* ----------- BEGIN: Code under measurement --------- */
  start_timer(t);
  int is_valid = user_verify(&usk, &ipk);
  time = stop_timer(t);
  /* ----------- END: Code under measurement ----------- */
  if (!is_valid) {
    printf("FATAL ERROR: benchmarked user key is not valid\n");
  }
  usk_clear(&usk);
  issuer_keys_clear(&ipk, &isk);
  return time;
}

double user_verify_invalid_bench(timer* t) {
  double time;
  isk_t isk;
  ipk_t ipk;
  usk_t usk;
  uint8_t state[STATE_BYTES];

  issuer_keys_init(&ipk, &isk);
  usk_init(&usk);
  issuer_keygen(&ipk, &isk);
  randombytes(state, STATE_BYTES);
  issuer_sign(&usk, state, &isk, &ipk);
  poly_q_set_coeff(usk.v2[0], 0, (coeff_q)1 + poly_q_get_coeff_centered(usk.v2[0], 0));
  /* ----------- BEGIN: Code under measurement --------- */
  start_timer(t);
  int is_valid = user_verify(&usk, &ipk);
  time = stop_timer(t);
  /* ----------- END: Code under measurement ----------- */
  if (is_valid) {
    printf("FATAL ERROR: benchmarked tampered user key is valid\n");
  }
  usk_clear(&usk);
  issuer_keys_clear(&ipk, &isk);
  return time;
}