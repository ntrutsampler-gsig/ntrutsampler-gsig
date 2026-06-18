#ifndef ARITH_H
#define ARITH_H

/**********************************************
* This header serves as an umbrella header for all arithmetic backends.
**********************************************/
#include "poly_real.h"
#include "poly_q.h"
#include "random.h"

void arith_setup(void);
void arith_teardown(void);

/**********************************************
* All functions (such as conversions) which rely on multiple "arithmetics" go here.
* Function that work on different types (i.e. poly vs d-vector vs k-vector vs matrix) should
* be placed in the corresponding "arithmetics".
**********************************************/

void poly_real_from_poly_q(poly_real res, const poly_q arg);

void poly_q_from_poly_real(poly_q res, const poly_real arg);

void poly_q_samplefz(poly_q res, const poly_real f, const poly_real c);

#endif /* ARITH_H */
