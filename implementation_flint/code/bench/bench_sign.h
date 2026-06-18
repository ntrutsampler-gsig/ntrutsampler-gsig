#ifndef BENCH_SIGN_H
#define BENCH_SIGN_H

#include "benchmark.h"

double issuer_keygen_bench(timer* t);
double issuer_sign_bench(timer* t);
double user_verify_valid_bench(timer* t);
double user_verify_invalid_bench(timer* t);

#endif /* BENCH_SIGN_H */

