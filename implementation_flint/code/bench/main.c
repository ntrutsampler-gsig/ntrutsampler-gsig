#include <stdio.h>
#include <time.h>

#include "sign.h"
#include "randombytes.h"
#include "random.h"

#include "bench_sign.h"
#include "bench_sampler.h"

#define BENCH_ITERATIONS 10

int main(void) {
  arith_setup();
  random_init();

  printf("[+] Running ./build/bench with N = %d iterations\n\n", BENCH_ITERATIONS);
  printf("______________________________________________________________________________________________________________________________________________________________\n");

  printf("%-30s\t%26s%4s%27s\t%27s%6s%28s\n", "", "", "Time", "", "", "Cycles", "");
  printf("%-30s\t%57s\t%61s\n", "", "---------------------------------------------------------", "-------------------------------------------------------------");
  printf("%-30s\t%9s\t%9s\t%9s\t%9s\t%13s\t%13s\t%13s\t%13s\n", "Benchmarked Functionality", "mean (ms)", "med (ms)", "min (ms)", "max (ms)", "mean (cycles)", "med (cycles)", "min (cycles)", "max (cycles)");
  printf("______________________________________________________________________________________________________________________________________________________________\n");

  printf("\n");

  benchmark("ISSUER_KEYGEN", BENCH_ITERATIONS, issuer_keygen_bench);
  benchmark("ISSUER_SIGN", BENCH_ITERATIONS, issuer_sign_bench);
  benchmark("USER_VERIFY (Valid)", BENCH_ITERATIONS, user_verify_valid_bench);
  benchmark("USER_VERIFY (Invalid)", BENCH_ITERATIONS, user_verify_invalid_bench);
  printf("\n");
  benchmark("PERTURBATION_SAMPLER", BENCH_ITERATIONS, perturbation_sampler_bench);
  benchmark("GADGET_SAMPLER", BENCH_ITERATIONS, gadget_sampler_bench);
  benchmark("NTRU_TSAMPLER", BENCH_ITERATIONS, ntru_tsampler_bench);

  printf("______________________________________________________________________________________________________________________________________________________________\n");


  arith_teardown();
  return 0;
}


