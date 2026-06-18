#ifndef POLY_REAL_H
#define POLY_REAL_H

#include <assert.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "params.h"
#include "fpr.h"

typedef fpr poly_real[PARAM_N];

typedef fpr coeff_real;

void poly_real_init(poly_real res);
void poly_real_clear(poly_real res);
void poly_real_zero(poly_real res);
void poly_real_set(poly_real res, const poly_real arg);
void poly_real_mul_scalar(poly_real res, const poly_real arg, const coeff_real f);
void poly_real_add_constant(poly_real res, const coeff_real c);
void poly_real_invert(poly_real res, const poly_real arg);
void poly_real_mul(poly_real res, const poly_real lhs, const poly_real rhs);
void poly_real_add(poly_real res, const poly_real lhs, const poly_real rhs);
void poly_real_sub(poly_real res, const poly_real lhs, const poly_real rhs);
void poly_real_set_si(poly_real res, size_t n, int64_t c);
void poly_real_samplefz(poly_real res, const poly_real f, const poly_real c);
int64_t poly_real_get_coeff_rounded(const poly_real arg, size_t n);

#endif /* POLY_REAL_H */
