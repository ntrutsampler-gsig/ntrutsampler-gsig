#ifndef RANDOM_H
#define RANDOM_H

#include <inttypes.h>
#include <x86intrin.h>
#include "fpr.h"

/*
	Code from random_aesni.c
*/

//public API
void random_init(void);

int sampleZ_small_wdth(fpr mu, fpr isigma);
int64_t SampleZ(fpr c, fpr sigma);

#endif /* RANDOM_H */
