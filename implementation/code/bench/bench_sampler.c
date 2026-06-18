#include "bench_sampler.h"

#include "sign.h"
#include "randombytes.h"
#include "random.h"
#include "arith.h"
#include "sampling.h"
#include "poly_q_sampling.h"

double perturbation_sampler_bench(timer* t) {
  double time;
  isk_t isk;
  ipk_t ipk;
  poly_q p1, p2[PARAM_KH];

  issuer_keys_init(&ipk, &isk);
  poly_q_init(p1);
  for (size_t i = 0; i < PARAM_KH; i++) {
    poly_q_init(p2[i]);
  }

  issuer_keygen(&ipk, &isk);
  /* ----------- BEGIN: Code under measurement --------- */
  start_timer(t);
  perturbation_sampler(p1, p2, isk.g, isk.f_prime, isk.s_p);
  time = stop_timer(t);
  /* ----------- END: Code under measurement ----------- */
  issuer_keys_clear(&ipk, &isk);
  poly_q_clear(p1);
  for (size_t i = 0; i < PARAM_KH; i++) {
    poly_q_clear(p2[i]);
  }
  return time;
}

double gadget_sampler_bench(timer* t) {
  double time;
  size_t i;
  isk_t isk;
  ipk_t ipk;
  poly_q tag, taginv;
  poly_q w, zL, zH[PARAM_KH];
  uint8_t state[STATE_BYTES];

  poly_q_init(w);
  poly_q_init(zL);
  poly_q_init(tag);
  poly_q_init(taginv);
  for (i = 0; i < PARAM_KH; i++) {
    poly_q_init(zH[i]);
  }
  issuer_keys_init(&ipk, &isk);

  issuer_keygen(&ipk, &isk);
  randombytes(state, STATE_BYTES);
  poly_q_uniform(w, ipk.seed, DOMAIN_SEPARATOR_U);
  tag_gen(tag, state);
  poly_q_invert_mod_bH(taginv, tag);
  /* ----------- BEGIN: Code under measurement --------- */
  start_timer(t);
  gadget_sampler(zL,zH,w,tag,taginv,isk.f,isk.finv_bH);
  time = stop_timer(t);
  /* ----------- END: Code under measurement ----------- */
  poly_q_clear(w);
  poly_q_clear(zL);
  poly_q_clear(tag);
  poly_q_clear(taginv);
  for (i = 0; i < PARAM_KH; i++) {
    poly_q_clear(zH[i]);
  }
  issuer_keys_clear(&ipk, &isk);

  return time;
}

double ntru_tsampler_bench(timer* t) {
  double time;
  size_t i;
  isk_t isk;
  ipk_t ipk;
  poly_q tag, taginv, u, v1, v2[PARAM_KH];
  uint8_t state[STATE_BYTES];

  issuer_keys_init(&ipk, &isk);
  poly_q_init(tag);
  poly_q_init(taginv);
  poly_q_init(u);
  poly_q_init(v1);
  for (i = 0; i < PARAM_KH; i++) {
    poly_q_init(v2[i]);
  }
 
  issuer_keygen(&ipk, &isk);
  poly_q_uniform(u, ipk.seed, DOMAIN_SEPARATOR_U);
  randombytes(state, STATE_BYTES);
  tag_gen(tag, state);
  poly_q_invert_mod_bH(taginv, tag);
  /* ----------- BEGIN: Code under measurement --------- */
  start_timer(t);
  ntru_tsampler(v1, v2, ipk.h, u, tag, taginv, &isk);
  time = stop_timer(t);
  /* ----------- END: Code under measurement ----------- */
  issuer_keys_clear(&ipk, &isk);
  poly_q_clear(tag);
  poly_q_clear(taginv);
  poly_q_clear(u);
  poly_q_clear(v1);
  for (i = 0; i < PARAM_KH; i++) {
    poly_q_clear(v2[i]);
  }

  return time;
}