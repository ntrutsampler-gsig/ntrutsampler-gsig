#ifndef PARAMS_H
#define PARAMS_H

/*************************************************
* Domain separators for XOF expansion
**************************************************/
#define DOMAIN_SEPARATOR_F 0
#define DOMAIN_SEPARATOR_E 1
#define DOMAIN_SEPARATOR_A3 2
#define DOMAIN_SEPARATOR_U 3

// Length of the public and secret seeds
#define SEED_BYTES 32
// Length of the state
#define STATE_BYTES 64

/*************************************************
* NTRU-TSampler parameters
**************************************************/
// Ring degree for the sampler
#define PARAM_N 1024
// Log of ring degree for the sampler
#define PARAM_LOG_N 10
// Ring degree for the proof
#define PARAM_N_PI 64
// Ring degree gap between R and R_pi
#define PARAM_K_PI 16
// Modulus for the sampler
#define PARAM_Q 123877L
// Low modulus
#define PARAM_QL 733L
// High modulus
#define PARAM_QH 169L
// Modulus inverse
#define PARAM_QL_MUL_QL_INVMOD_QH 63038L
#define PARAM_QH_MUL_QH_INVMOD_QL 60840L
// Modulus bit-length upper bound for uniform sampling
#define PARAM_Q_BITLEN 17
// G_H dimension
#define PARAM_KH 2
// G_H base
#define PARAM_BH 13
// Hamming weight of the tags
#define PARAM_W 8
// Squared verification bound on v_1
#define PARAM_B1SQ 5980214451L
// Squared verification bound on [v_2 | v_3]
#define PARAM_B2SQ 8909931057L

// Bound on square spectral norm of e^T
#define PARAM_BE_POW_4 35620716.80454716086387634277
// Bound on spectral norm of the f_i
#define PARAM_BFSQ 4096.00000000000000000000
// Difference s_1² - s_L²
#define PARAM_S1SQ_SLSQ 23283790.35054315999150276184
// Squared Gaussian parameter s_2^2
#define PARAM_S2SQ 15979459.46655039116740226746
// Gaussian parameter s_2
#define PARAM_S2 3997.43160873959050150006
// Negated squared Gaussian parameter s_H^2
#define PARAM_NEG_SHSQ -1938.94565414150201831944
// Gaussian parameter s_L/q_L for Z-sampling of zL
#define PARAM_SL_DIV_QL 3.35446204446931162124
// Gaussian parameter s_H/b_H for Z-sampling of zH
#define PARAM_SH_DIV_BH 3.38718927843394679300

// Negated inverse of q_L : -1/q_L
#define PARAM_NEG_INV_QL -0.00136425648021828104
// Negated inverse of b_H : -1/b_H
#define PARAM_NEG_INV_BH -0.07692307692307692735

// 0
#define PARAM_ZERO 0
// 1
#define PARAM_ONE 1.0
// 1/2
#define PARAM_ONE_HALF 0.5
// 6
#define PARAM_SIX 6.0
// UINT16_MAX (65535)
#define PARAM_UINT16_MAX 65535.0
// Inverse of sqrt(2*pi)
#define PARAM_INV_SQRT_2PI 0.39894228040143270286


/*************************************************
* Small widths sampler constants
**************************************************/
/*
// Inverse of s_L/q_L scaled by sqrt(2*.pi) : 1/((s_L/q_L)/sqrt(2.pi)) = sqrt(2*pi)*q_L / s_L
#define PARAM_INV_SL_DIV_QL 0.74725194126546101714
// Inverse of s_H/b_H scaled by sqrt(2*.pi) : 1/((s_H/b_H)/sqrt(2.pi)) = sqrt(2*pi)*b_H / s_H
#define PARAM_INV_SH_DIV_BH 0.74003194642548286719
// Inverse of log(2)
#define PARAM_INV_LOG2 1.44269504088896338700
// Natural log of 2
#define PARAM_LOG2 0.69314718055994528623
// Two power 63
#define PARAM_2_POW_63 9223372036854775808L
// Inverse of two power 63
#define PARAM_INV_2_POW_63 -0.00000000000000000011
// Inverse of 2*sigma_max^2
#define PARAM_INV_2_SIGMAMAX_SQ 0.15086504887537272035
// Sigma_min for small width sampler
#define PARAM_SIGMAMIN 1.29828033434429190862
*/

#endif /* PARAMS_H */
