#ifndef BENCH_SAMPLER_H
#define BENCH_SAMPLER_H

#include "benchmark.h"

double perturbation_sampler_bench(timer* t);
double gadget_sampler_bench(timer* t);
double ntru_tsampler_bench(timer* t);

#endif /* BENCH_SAMPLER_H */

