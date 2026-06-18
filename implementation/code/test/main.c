#include <stdio.h>
#include "poly_q_sampling.h"
#include "arith.h"
#include "sign.h"
#include "sampling.h"
#include "precomputations.h"
#include "randombytes.h"
#include "random.h"

#define NTESTS 1
#define NSUBTESTS 5

static int issuer_keygen_test(void) {
  size_t i,j;
  int rval = 1;
  isk_t isk;
  ipk_t ipk;
  poly_q tmp[PARAM_KH];

  printf("\nissuer_keygen_test\n");

  issuer_keys_init(&ipk, &isk);
  for (i = 0; i < PARAM_KH; i++) {
    poly_q_init(tmp[i]);
  }

  for (i = 0; i < NSUBTESTS; i++) {
    issuer_keygen(&ipk, &isk);
    for (j = 0; j < PARAM_KH; j++) {
      if (f_sq_spectral_norm(isk.f[j])) {
        printf("issuer_keygen generated keys with f_j with large spectral norm.\n");
        rval = 0;
        goto issuer_keygen_test_cleanup;
      }
    }
    if (e_sq_spectral_norm(tmp, isk.e)) {
      printf("issuer_keygen generated keys with e*.e with large spectral norm.\n");
      rval = 0;
      goto issuer_keygen_test_cleanup;
    }
    
    for (j = 0; j < PARAM_KH; j++) {
      poly_q_invert_mod_q_bH(tmp[0], tmp[1], isk.f[j]);
      poly_q_mul(tmp[0], isk.e[j], tmp[0]);
      if (!poly_q_equal(tmp[0], ipk.h[j])) {
        printf("issuer_keygen generated keys with h != e.f^{-1} mod q.\n");
        rval = 0;
        goto issuer_keygen_test_cleanup;
      }
    }

    poly_q_set_coeff(isk.e[0], PARAM_N/2, 3); // changing one coefficient of e
    poly_q_invert_mod_q_bH(tmp[0], tmp[1], isk.f[0]);
    poly_q_mul(tmp[0], isk.e[0], tmp[0]);
    if (poly_q_equal(tmp[0], ipk.h[0])) {
      printf("found e' such that h = e'.f^{-1} mod q.\n");
      rval = 0;
      goto issuer_keygen_test_cleanup;
    }

    for (j = 0; j < PARAM_KH; j++) {
      poly_q_mul(tmp[0], isk.f[j], isk.finv_bH[j]);
      poly_q_mod_bH(tmp[1], tmp[0]);
      if (!nmod_poly_is_one(tmp[1])) {
        printf("inverse modulo bH is not correct.\n");
        rval = 0;
        goto issuer_keygen_test_cleanup;
      }
    }

    printf(":");
    fflush(stdout);
  }

issuer_keygen_test_cleanup:
  issuer_keys_clear(&ipk, &isk);
  for (i = 0; i < PARAM_KH; i++) {
    poly_q_clear(tmp[i]);
  }

  return rval;
}

static int gadget_sampler_test(void) {
  size_t i,j,k;
  int rval = 1;
  coeff_q bexpi;
  isk_t isk;
  ipk_t ipk;
  poly_q tag, taginv, tmp, w, zL, zH[PARAM_KH];
  uint8_t state[STATE_BYTES];
  
  printf("\ngadget_sampler_test\n");

  poly_q_init(w);
  poly_q_init(tmp);
  poly_q_init(zL);
  poly_q_init(tag);
  poly_q_init(taginv);
  for (i = 0; i < PARAM_KH; i++) {
    poly_q_init(zH[i]);
  }
  issuer_keys_init(&ipk, &isk);

  for (i = 0; i < NSUBTESTS; i++) {
    issuer_keygen(&ipk, &isk);
    randombytes(state, STATE_BYTES);
    poly_q_uniform(w, ipk.seed, DOMAIN_SEPARATOR_U);
    for (j = 0; j < NSUBTESTS; j++) {
      tag_gen(tag, state);
      poly_q_invert_mod_bH(taginv, tag);

      gadget_sampler(zL,zH,w,tag,taginv,isk.f,isk.finv_bH);

      bexpi = (coeff_q)PARAM_QL;
      for (k = 0; k < PARAM_KH; k++) {
        poly_q_mul_scalar(tmp, tag, bexpi);
        poly_q_mul(tmp, tmp, isk.f[k]);
        poly_q_mul(tmp, tmp, zH[k]);
        poly_q_add(zL, zL, tmp);
        bexpi *= PARAM_BH;
      }

      if (!poly_q_equal(zL, w))
      {
        printf("gadget_sampler returned [zL | zH] that do not verify the correct equation.\n");
        rval = 0;
        goto gadget_sampler_test_cleanup;
      }
      printf(":");
      fflush(stdout);
    }
  }

gadget_sampler_test_cleanup:
  poly_q_clear(w);
  poly_q_clear(tmp);
  poly_q_clear(zL);
  poly_q_clear(tag);
  poly_q_clear(taginv);
  for (i = 0; i < PARAM_KH; i++) {
    poly_q_clear(zH[i]);
  }
  issuer_keys_clear(&ipk, &isk);
  return rval;
}

static int user_key_test(void) {
  size_t i,j;
  int rval = 1;
  isk_t isk;
  ipk_t ipk;
  usk_t usk;
  uint8_t state[STATE_BYTES];
  
  printf("\nuser_key_test\n");

  issuer_keys_init(&ipk, &isk);
  usk_init(&usk);

  for (i = 0; i < NSUBTESTS; i++) {
    issuer_keygen(&ipk, &isk);
    randombytes(state, STATE_BYTES);
    for (j = 0; j < NSUBTESTS; j++) {
      issuer_sign(&usk, state, &isk, &ipk);
      if (!user_verify(&usk, &ipk))
      {
        printf("user_verify returned zero for a valid user key.\n");
        rval = 0;
        goto user_key_test_cleanup;
      }

      poly_q_set_coeff(usk.v2[0], 0, (coeff_q)1 + poly_q_get_coeff_centered(usk.v2[0], 0));
      if (user_verify(&usk, &ipk))
      {
        printf("user_verify returned non-zero for a tampered user key.\n");
        rval = 0;
        goto user_key_test_cleanup;
      }
      printf(":");
      fflush(stdout);
    }
  }

user_key_test_cleanup:
  usk_clear(&usk);
  issuer_keys_clear(&ipk, &isk);
  return rval;
}

int main(void) {
  int pass = 1;
  arith_setup();
  random_init();
  printf("Hello from the unit tests.\n");
  for (int i = 0; i < NTESTS; i++)
  {
    pass &= issuer_keygen_test();
    pass &= gadget_sampler_test();
    pass &= user_key_test();

    if (!pass)
    {
      printf("FAILED!\n");
      break;
    } else {
      printf(".");
    }
  }
  if (pass)
  {
    printf("\npassed.\n");
  }
  arith_teardown();
  return 0;
}
